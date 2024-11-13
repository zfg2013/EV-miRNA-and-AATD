rm(list=ls())

library(tidyverse)
library(edgeR)
library(pROC)
library(readxl)


WD_root <- "/Users/zhuo/Dropbox/research/Khodayari/AATD_miRNA"
WD_bioinformatics <- file.path(WD_root, "bioinformatics")
WD_data <- file.path(WD_bioinformatics, "countMatrix", "bowtie2_secondary")
WD_data_batch2 <- file.path(WD_bioinformatics, "countMatrix", "bowtie2_secondary_batch2")
WD_phenotype_batch2 <- file.path(WD_bioinformatics, "countMatrix", "phenotype_batch2")
dir(WD_phenotype_batch2)

WD_statAnalysis <- file.path(WD_root, "statAnalysis_isotype")
WD_DE <- file.path(WD_statAnalysis, "DE_PASD_01_23")
WD_roc <- file.path(WD_DE, "roc")

dir.create(WD_roc)

afile <- paste0("DE_","Table", ".csv")
DEresult2 <- read_csv(file.path(WD_DE, afile))
pheno_edata <- read_csv(file.path(WD_DE, "phenotype_inuse.csv"))
pheno_batch2 <- read_excel(file.path(WD_phenotype_batch2, "summary_bowtie2_batch2_withID.xlsx"))
pheno_batch2 <- pheno_batch2 %>% 
	separate(sample, into=c("ID", "i7"), sep="_") %>%
	mutate(PASD_binary=ifelse(PASD>=2, 1, 0))

data_count0 <- get(load(file.path(WD_data, "count_no.rdata")))
rowSums <- apply(data_count0,1,function(x) sum(x))
data_count_expressed <- data_count0[rowSums > 0, ]
adata <- data_count_expressed[, pheno_edata$ID]
data_cpm <- cpm(adata)
data_cpm_log2 <- log2(data_cpm + 1)


data_count0_batch2 <- get(load(file.path(WD_data_batch2, "count_no.rdata")))
data_cpm_batch2 <- cpm(data_count0_batch2)
data_cpm_log2_batch2 <- log2(data_cpm_batch2 + 1)

stopifnot(pheno_batch2$ID == colnames(data_cpm_log2_batch2))


group <- pheno_edata$PASD_binary
type <- ifelse(group==0, "01", "23")

group_batch2 <- ifelse(pheno_batch2$PASD>=2, 1, 0)
type_batch2 <- ifelse(group_batch2==0, "01", "23")

 miR-23a-3p, miR-223-3p, miR-30a-3p, miR-374a-3p, and miR-191-5p
miRNA_panels <- c(
    "hsa-miR-23a-3p",
    "hsa-miR-223-3p",
    "hsa-miR-30a-5p",
    "hsa-miR-374a-3p",
    "hsa-miR-191-5p"
    )

len_miRNAs <- length(miRNA_panels)
combn_index <- combn(len_miRNAs,4)

miRNA7_all <- unique(unlist(miRNA_panels))

miRNA7_all %in% rownames(data_cpm_log2) 

panels_batch1 <- t(scale(t(data_cpm_log2[miRNA7_all,])))
panels_batch2 <- t(scale(t(data_cpm_log2_batch2[miRNA7_all,])))

WD_otherPanel <- file.path(WD_roc, "otherPanel_standardize_v4")
dir.create(WD_otherPanel,re=TRUE)
setwd(WD_otherPanel)

write.csv(panels_batch1, "data_panels_batch1.csv")
write.csv(panels_batch2, "data_panels_batch2.csv")


i=1

combn_index <- combn(len_miRNAs,4)

ROC_result_combine <- NULL


for(i in 1:ncol(combn_index)){
	aname <- i
	miRNA_panel <- miRNA_panels[combn_index[,i]]

  adata_7miRNA <-  as_tibble(t(data_cpm_log2[miRNA_panel,]))
  adata_7miRNA2 <- tibble(ID=colnames(data_cpm_log2), adata_7miRNA)
  pheno_edata_miRNA <- pheno_edata %>% 
  	left_join(adata_7miRNA2) 
  pheno_edata_miRNA2 <- pheno_edata_miRNA %>% 
  	select(PASD_binary, starts_with("hsa"))
  
  ##
  alogistic <- glm(PASD_binary ~ ., family = binomial(link = "logit"), data = pheno_edata_miRNA2)
  pheno_edata_miRNA3 <- pheno_edata_miRNA %>% mutate(fitted = predict(alogistic, pheno_edata_miRNA2, type="response"))
	pheno_edata_miRNA3$fitted

  aroc <- with(pheno_edata_miRNA3, roc(PASD_binary, fitted))
  anindex <- which.max(aroc$sens * aroc$spec)
  asens <- aroc$sens[anindex]
  aspec <- aroc$spec[anindex]
  aauc <- aroc$auc
  asens_spec <- asens * aspec
	athreshold <- aroc$thresholds[anindex]
	if(aroc$direction == "<"){
		case <- 1
	} else {
		case <- 0
	}
	
  aROC_result_batch1 <- tibble(miRNA1 = miRNA_panel[1],
      miRNA2 = miRNA_panel[2],
  	miRNA3 = miRNA_panel[3],
	miRNA4 = miRNA_panel[4],
	sensitivity_batch1 = asens,
	specificity_batch1 = aspec, 
	AUC_batch1 = aauc,
	`sensitivity*specificity_batch1` = asens_spec)
		
	
	## batch2
  adata_7miRNA_batch2 <-  as_tibble(t(data_cpm_log2_batch2[miRNA_panel,]))
  adata_7miRNA2_batch2 <- tibble(ID=colnames(data_cpm_log2_batch2), adata_7miRNA_batch2)
  pheno_edata_miRNA_batch2 <- pheno_batch2 %>% 
  	left_join(adata_7miRNA2_batch2) 
  pheno_edata_miRNA2_batch2 <- pheno_edata_miRNA_batch2 %>% 
  	select(PASD_binary, starts_with("hsa"))

	#
  pheno_edata_miRNA3_batch2 <- pheno_edata_miRNA2_batch2 %>% mutate(fitted = predict(alogistic, pheno_edata_miRNA2_batch2, type="response"))
	pheno_edata_miRNA3$fitted
	pheno_edata_miRNA3_batch2$fitted
	
	
  aroc <- with(pheno_edata_miRNA3_batch2, roc(PASD_binary, fitted))
	alogic <- factor(ifelse(pheno_edata_miRNA3_batch2$fitted > athreshold, case, 1-case), level=c(0,1))
  atable <- table(alogic, pheno_edata_miRNA3_batch2$PASD_binary)
	asens <- ifelse(atable[2,2]==0, 0, atable[2,2]/sum(atable[,2]))
	aspec <- ifelse(atable[1,1]==0, 0, atable[1,1]/sum(atable[,1]))			
  aauc <- aroc$auc
  asens_spec <- asens * aspec

  aROC_result_batch2 <- tibble(miRNA1 = miRNA_panel[1],
															miRNA2 = miRNA_panel[2],
															miRNA3 = miRNA_panel[3],
															sensitivity_batch2 = asens,
															specificity_batch2 = aspec, 
															AUC_batch2 = aauc,
															`sensitivity*specificity_batch2` = asens_spec)
	
	#
	aROC_result_batch <- aROC_result_batch1 %>% left_join(aROC_result_batch2)
	
	ROC_result_combine <- rbind(ROC_result_combine, aROC_result_batch)
	
	
}


write.csv(ROC_result_combine, "ROC_result_batch_cobmine4.csv")
