rm(list=ls())

library(tidyverse)
library(edgeR)
library(pROC)

WD_root <- "/Users/zhuo/Dropbox/research/Khodayari/AATD_miRNA"
WD_bioinformatics <- file.path(WD_root, "bioinformatics")
WD_data <- file.path(WD_bioinformatics, "countMatrix", "bowtie2_secondary")

WD_statAnalysis <- file.path(WD_root, "statAnalysis_isotype")
WD_DE <- file.path(WD_statAnalysis, "DE_PASD_01_23")
WD_roc <- file.path(WD_DE, "roc")

dir.create(WD_roc)

afile <- paste0("DE_","Table", ".csv")
DEresult2 <- read_csv(file.path(WD_DE, afile))
pheno_edata <- read_csv(file.path(WD_DE, "phenotype_inuse.csv"))

data_count0 <- get(load(file.path(WD_data, "count_no.rdata")))
rowSums <- apply(data_count0,1,function(x) sum(x))
data_count_expressed <- data_count0[rowSums > 0, ]

adata <- data_count_expressed[, pheno_edata$ID]
data_cpm <- cpm(adata)
data_cpm_log2 <- log2(data_cpm + 1)

group <- pheno_edata$PASD_binary
type <- ifelse(group==0, "01", "23")

DEresult3 <- DEresult2 %>% 
    mutate(sensitivity = NA, specificity = NA, AUC = NA, `sensitivity*specificity` = NA)
for(i in 1:nrow(DEresult2)){
    amiRNA <- DEresult2$gene_id[i]
    agene <- data_cpm_log2[amiRNA,]
    aroc <- roc(type, agene)
    anindex <- which.max(aroc$sens * aroc$spec)
    asens <- aroc$sens[anindex]
    aspec <- aroc$spec[anindex]
    aauc <- aroc$auc
    asens_spec <- asens * aspec
    
    DEresult3$sensitivity[i] <- asens
    DEresult3$specificity[i] <- aspec
    DEresult3$AUC[i] <- aauc
    DEresult3$`sensitivity*specificity`[i] <- asens_spec
    
}

DEresult4 <- DEresult3 %>% arrange(desc(`sensitivity*specificity`))
## visualization
setwd(WD_roc)
write.csv(DEresult4, "DE_roc.csv", row.names = FALSE)


topN <- 39
for(i in 1:topN){
	amiRNA <- DEresult4$gene_id[i]
  agene <- data_cpm_log2[amiRNA,]
  aroc <- roc(type, agene)

	atitle <- paste0(amiRNA, ". AUC = ", signif(aroc$auc, 2))
	afile <- paste0("ROC_top", i, ".pdf")
	pdf(afile,5,5)	
	plot(x = 1 - aroc$spec, y = aroc$sens, pch = 4, col = "blue", type = "n",
	     xlab = "1 - Specificity", ylab = "Sensitivity", cex = 2, 
		 xlim = c(0, 1), ylim = c(0, 1), main = atitle)
	lines(x = 1 - aroc$spec, y = aroc$sens, type = "l", col = "blue",
	     xlab = "1 - Specificity", ylab = "Sensitivity", cex = 2, 
		 xlim = c(0, 1), ylim = c(0, 1))
	#
	dev.off()

}
