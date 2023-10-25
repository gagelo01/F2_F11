#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(openxlsx)
wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/F2_finngen"
setwd(wd)
res_unicis <- fread ( "Data/Modified/res_unicis.txt")
res_sign <- fread( "Data/Modified/res_sign.txt")
res_unicis[clean_variable_name=="", clean_variable_name:=outcome]
res_multicis <- fread("Data/Modified/res_multicis_meta.txt")
dt_pan <-fread("Data/Modified/dt_pan.txt")
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
efficacy_outcome <- c("dis-1-1", "dis-14-7", "dis-14-8","dis-14-9", "dis-14-10")
safety_outcome <- c("dis-15-157", "dis-15-1516", "dis-15-2495", "dis-15-2844", "trait-7-1")
####
dataset <- df_index[id %in% c(efficacy_outcome, safety_outcome), ]
res_multicis[,note:=NULL]
res_multicis[,lci:=b-se*1.96]
res_multicis[,uci:=b+se*1.96]
res_pan <- merge(dt_pan, dt_gene_region, by.x = "id.exposure", by.y = "id")

list_supdat <- list( "Supplementary Table 1" = dataset,
                     "Supplementary Table 2" = res_unicis,
                     "Supplementary Table 3" = res_multicis,
                     "Supplementary Table 4" = res_pan)
toremove<-c("cochranQ", "cochranQ_p", "gene_region", "vcffile")
list_supdat[2:3] <- map(list_supdat[2:3], function(x) {
  x <- x[id.outcome %in% c(efficacy_outcome, safety_outcome), ]
  x[, (toremove) := NULL]
})

#

dt_title <- data.table(title = paste0("ST", c(1:4)),
                       caption = c( "Description of the datasets used for Two-Sample Mendelian randomization.",
                                    "Results of the uni-cis MR investigation for the effect blood protein levels and hepatic gene levels on efficacy and safety outcomes. Colocalisation results are also presented.",
                                    "Results of the multi-cis MR investigation for the effect blood protein levels on efficacy and safety outcomes.",
                                     "Results of the pan MR investigation for the effect blood protein levels on efficacy and safety outcomes."))




###
col_description<- vector(mode = "list", length = length(list_supdat))
legends = c("trait = Unique trait identifier; group_name = the name of the study group;  year = year data was published; author = author of the data; consortium = consortium for the sample; sex = sex of the sample (in this study always Males and Females); population = Ancestry of the sample (in this study always European); unit = The unit either standard deviation (SD) or log Odds ratio (log(orr)); nsnp = number of SNPs in the summary statistics; sample_size = The maximum sample size; ncase = The maximum number of cases; ncontrol = the maximum number of controls; pmid = Pubmed ID",
            "outcome = unique name of the outcome; exposure = unique name of the exposure; method = name of the method to estimate the causal effect of the exposure on the outcome; nsnp = The number of SNPs genetic instruments used to compute the estimate or the name of the SNP if only one; b = the estimate scaled per SD or log(orr);  se =  standard error of the estimate; pval = the p-value of the estimate; lci = the 95% confidence intervall lower bound of the effect; uci = the 95% confidence intervall upper bound of the effect; fstat the F-statistic;rsq = the variance explaine; steiger_dir = the steiger test; steiger_pval = the steiger test pval, posprob_coloc_PPH0:4 = coloc posterior probabilitoy of hypothesis 0 to 4; posprob_colocH4.SNP = the SNP prioritised to be causal for both traits;hgnc = the hgnc gene name; UniProt the UniProt protein name; study = The study from which the exposure is taken; clean_variable_name = The outcome name cleaner; outcome_category = The outcome category; is_altering_variant = Does the variant alter the amino acids sequence of the protein ? Annotated with Variant Effect Predictor",                
            "outcome = unique name of the outcome; exposure = unique name of the exposure; method = name of the method to estimate the causal effect of the exposure on the outcome; nsnp = The number of SNPs genetic instruments used to compute the estimate; b = the estimate scaled per SD or log(orr);  se =  standard error of the estimate; pval = the p-value of the estimate; lci = the 95% confidence intervall lower bound of the effect; uci = the 95% confidence intervall upper bound of the effect; study = The study from which the exposure is taken; clean_variable_name = The outcome name cleaner; outcome_category = The outcome category; cochranQ.multi_cis = the value for the cochran'Q; cochranQpval.multi_cis = the pvalue of the cochranQ",                
            "outcome = unique name of the outcome; exposure = unique name of the exposure; method = name of the method to estimate the causal effect of the exposure on the outcome; nsnp = The number of SNPs genetic instruments used to compute the estimate; b = the estimate scaled per SD or log(orr);  se =  standard error of the estimate; pval = the p-value of the estimate; lci = the 95% confidence intervall lower bound of the effect; uci = the 95% confidence intervall upper bound of the effect; type_of_test = categorisation of the method based on Slob and Burgess pmid = 32249995.study = The study from which the exposure is taken; clean_variable_name = The outcome name cleaner; outcome_category = The outcome category; hgnc = the hgnc gene name; UniProt the UniProt protein name; study = The study from which the exposure is taken;")
            
col_description[[1]] <- tribble(
  ~x, ~y,  
  "id", "a unique id", 
  "trait", "Unique trait identifier",
  "year", "the year the GWAS was published", 
  "trait", "The author of the GWAS",
  "consortium", "the name of the consortium of the GWAS",
  "sex", "sex included",
  "population", "ancestry",
  "sample_size", "the sample_size",
  "pmid", "the pubmed ID",
  "ncase", "the number of cases",
  "ncontrol", "the number of controls",
  "adjustments", "what variables were included in the GWAS adjustment model",
  "clean_variable_name", "A publication ready name that can be used to plot figures."
) %>% as.data.table(.)

col_description[[2]] <- tribble(
  ~x, ~y,  
  "id.outcome", "a unique id for the outcome", 
  "id.exposure", "a unique id for the exposure",
  "outcome", "a unique name for the outcome", 
  "exposure", "a unique name for twosample MR",
  "method", "name of the method to estimate the causal effect of the exposure on the outcome",
  "nsnp", "The number of SNPs genetic instruments used to compute the estimate or the name of the SNP if only one",
  "b", "the estimate scaled per SD or log(orr);  se =  standard error of the estimate",
  "se", "standard error",
  "pval", "p-value",
  "lci", "lower confidence interval",
  "uci", "upper confidence interval",
  "fstat", "he F-statistic",
  "rsq", "The variance explained",
  "steiger_dir", "the steiger directionnality test",
  "steiger_pval", "he steiger test pval",
  "posprob_coloc_PPH0:4",  "coloc posterior probabilitoy of hypothesis 0 to 4",
  "posprob_colocH4.SNP", "the SNP prioritised to be causal for both traits",
  "hgnc", "the hugo gene nomenclature gene name",
  "UniProt", "the UniProt protein name",
  "study", "The study from which the exposure is taken",
  "clean_variable_name", "A publication ready outcome name that can be used to plot figures.",
  "outcome_category", "The outcome category",
  "is_altering_variant", "Does the variant alter the amino acids sequence of the protein ? Annotated with Variant Effect Predictor", 
) %>% as.data.table(.)


col_description[[3]] <- tribble(
  ~x, ~y,  
  "id.outcome", "a unique id for the outcome", 
  "id.exposure", "a unique id for the exposure",
  "outcome", "a unique name for the outcome", 
  "exposure", "a unique name for twosample MR",
  "method", "name of the method to estimate the causal effect of the exposure on the outcome",
  "nsnp", "The number of SNPs genetic instruments used to compute the estimate or the name of the SNP if only one",
  "b", "the estimate scaled per SD or log(orr);  se =  standard error of the estimate",
  "se", "standard error",
  "pval", "p-value",
  "lci", "lower confidence interval",
  "uci", "upper confidence interval",
  "cochranQ.multi_cis", "the value for the cochran'Q",
  "cochranQpval.multi_cis", "the pvalue of the cochranQ",
  "hgnc", "the hugo gene nomenclature gene name",
  "UniProt", "the UniProt protein name",
  "study", "The study from which the exposure is taken",
  "clean_variable_name", "A publication ready outcome name that can be used to plot figures."
) %>% as.data.table(.)

col_description[[4]] <- tribble(
  ~x, ~y,  
  "id.outcome", "a unique id for the outcome", 
  "id.exposure", "a unique id for the exposure",
  "outcome", "a unique name for the outcome", 
  "exposure", "a unique name for twosample MR",
  "method", "name of the method to estimate the causal effect of the exposure on the outcome",
  "nsnp", "The number of SNPs genetic instruments used to compute the estimate or the name of the SNP if only one",
  "b", "the estimate scaled per SD or log(orr);  se =  standard error of the estimate",
  "se", "standard error",
  "pval", "p-value",
  "lci", "lower confidence interval",
  "uci", "upper confidence interval",
  "type_of_test", "categorisation of the method based on Slob and Burgess pmid = 32249995.study",
  "hgnc", "the hugo gene nomenclature gene name",
  "UniProt", "the UniProt protein name",
  "study", "The study from which the exposure is taken",
  "clean_variable_name", "A publication ready outcome name that can be used to plot figures.",
  "outcome_category", "The outcome category",
) %>% as.data.table(.)


bold_st <- createStyle(textDecoration = "Bold")
wb <- createWorkbook()
for(i in 1:length(list_supdat)) {
  addWorksheet(wb, sheetName =  dt_title[i, title])
  
  title <- dt_title[i,paste0(title, " : ", caption)]
  writeData(wb, sheet = i, x = title, startCol = 1, startRow = 1)
  addStyle(wb, sheet = i, style = bold_st, rows = 1, cols = 1:2)
  writeData(wb, sheet = i, x = col_description[[i]], startCol = 1, startRow = 2)
  addStyle(wb, sheet = i, style = bold_st, rows = 2:col_description[[i]][,.N+2], cols = 1)
  deleteData(wb, sheet = i, rows = 2, cols = 1:2, gridExpand = TRUE)
  writeData(wb, sheet = i, x = list_supdat[[i]], startCol = 1, startRow = col_description[[i]][,.N+4])
  addStyle(wb, sheet = i, style = bold_st, rows = col_description[[i]][,.N+4], cols = 1:length(colnames(list_supdat[[i]])), gridExpand = TRUE, stack = TRUE)
}
saveWorkbook(wb, file = "Results/supplementary_tables.xlsx", overwrite = TRUE)

message("This script finished without errors")

