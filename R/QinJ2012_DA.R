#Implement packages
library(pacman)
pacman::p_load("curatedMetagenomicData", "EnrichmentBrowser", "SummarizedExperiment", "dplyr", "stringr", "EnhancedVolcano", "dbplyr", "DESeq2", "png", "qusage")

#Format data for SummarizedExperiment

#Get meta
metahmp <- readr::read_tsv("C:/Users/Zachary Nachod/Documents/QinJ_2012_metadata.tsv")
hmp.metadata <- (metahmp)
hmp.metadata <- hmp.metadata[,-1]
hmp.metadata <- as.data.frame(hmp.metadata)
hmp.metadata["study_condition"][is.na(hmp.metadata["study_condition"])] <- 0

#reformat metadata
tempCategories<-colnames(hmp.metadata)
templist<-hmp.metadata$sample_id
hmp.metadata$sample_id<-NULL
colnames(hmp.metadata)<-NULL
hmp.metadata <- rbind(tempCategories[1:length(tempCategories)], hmp.metadata)
rownames(hmp.metadata)<-c("sample_id",templist)


hmp.metadata <- as.data.frame(hmp.metadata) %>%
  janitor::row_to_names(row_number = 1)

#Get bacterial abundance data
hmp<-curatedMetagenomicData("QinJ_20.+.relative_abundance", dryrun = FALSE) %>%
  mergeData()
assay<-hmp@assays@data@listData[["relative_abundance"]]
#assay<-t(assay)
#assay<-assay[2376:2671,]#index if multiple studies are combined
#assay<-t(assay)

#reformat bacteria names
reg_name <- function(MetaPhlAn){
  temp <- sub(".*g__", "", MetaPhlAn)
  gsub("\\|s__.*", "", temp)
  gsub("_noname", "", temp)
}

.getLast <- function(n)
{
  spl <- unlist(strsplit(n, "\\|"))
  spl[length(spl)]
}

rownames(assay) <- vapply(rownames(assay), reg_name,
                          character(1), USE.NAMES = FALSE)
rownames(assay) <- vapply(rownames(assay), .getLast,
                          character(1), USE.NAMES = FALSE)

mode(assay) <- "numeric"
hmp.se <- SummarizedExperiment(assays = list(assay), colData = hmp.metadata)

#get raw read counts
mode(hmp.se$number_reads) <- "numeric"
mode(assay(hmp.se)) <- "numeric"
hmp.counts = sweep(assay(hmp.se), 2, hmp.se$number_reads / 100, "*")
hmp.counts = round(hmp.counts)
hmp.se <- SummarizedExperiment(assays = hmp.counts, colData = hmp.metadata)

#location subset

#hmp.se <- subset(hmp.se, , hmp.se$location == "Boston")

grp1 = ifelse(hmp.se$study_condition == "T2D" , 1, 0)
hmp.se$GROUP = grp1

# keep only those microbes that have more than 430 read counts across all samples
keep <- rowSums(assay(hmp.se)) > 360 #change depending on size of
hmp.se <- hmp.se[keep,]

#differential expression analysis
assay(hmp.se) <- assay(hmp.se) + 1
hmp.se <- deAna(hmp.se, de.method = "DESeq2", filter.by.expr = F)

#reformat hmp microbe names, no MetaPhlAn string format
getGenus <- function(MetaPhlAn){
  temp <- gsub("s__", "", MetaPhlAn)
  gsub("_noname", "", temp)
  gsub("\\_.*","",temp)
}

rownames(hmp.se) <- vapply(rownames(hmp.se), getGenus,
                           character(1), USE.NAMES = FALSE)

#set-based enrichment analyses
bergeys.gmt <- qusage::read.gmt("C:/Users/Zachary Nachod/Documents/bug_physiologies.gmt")

#ORA
ora_results <- sbea("ora", hmp.se, bergeys.gmt, perm = 0)
ora_results <- gsRanking(ora_results, signif.only = F) %>%
  as.data.frame()
colnames(ora_results) <- sub("GENE", "MICROBE", colnames(ora_results))

#GSEA
gsea_results <- sbea("gsea", hmp.se, bergeys.gmt, perm = 1000)
gsea_results <- gsRanking(gsea_results, signif.only = F) %>%
  as.data.frame()
colnames(gsea_results) <- sub("GENE", "MICROBE", colnames(gsea_results))

gdata::keep(hmp.se, ora_results, gsea_results, sure = TRUE)

#Volcano plot for DESeq2 results
png(filename = "C:/Users/Zachary Nachod/Documents/data_raw_TEST/HMP_Volcano.png", res = 300, height = 1800, width = 1800)
EnhancedVolcano(rowData(hmp.se),
                lab = rownames(hmp.se),
                pCutoff = 0.05,
                title = 'Differential Abundance[DESeq2]',
                x = 'FC',
                y = 'ADJ.PVAL',
                xlim = c(-20, 20),
                pointSize = 1.5)
dev.off()

#Results output
write.csv(ora_results, file = "C:/Users/Zachary Nachod/Documents/data_raw_TEST/HMP_ora_res.csv")
write.csv(gsea_results, file = "C:/Users/Zachary Nachod/Documents/data_raw_TEST/HMP_gsea_res.csv")

#volcanoplot datatable
volcanoData <- as.data.frame(rowData(hmp.se))
volcanoDataGray <- filter(volcanoData, log2(abs(FC))<1, ADJ.PVAL > .05)
volcanoDataGreen <- filter(volcanoData, log2(abs(FC))>1, ADJ.PVAL > .05)
volcanoDataBlue <- filter(volcanoData, log2(abs(FC))<1, ADJ.PVAL < .05)
volcanoDataRed <- filter(volcanoData, log2(abs(FC))>1, ADJ.PVAL < .05)
