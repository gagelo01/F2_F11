#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(furrr)
library(metafor)

wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/F2_finngen"
setwd(wd)
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
gwasvcf::set_bcftools()
gwasvcf::set_plink()
options(future.globals.maxSize= 5e9)
plan(multicore, workers = 20, gc = TRUE)
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
dt_wald <- fread( "Data/Modified/dt_wald.txt")
dt_unicis <- fread( "Data/Modified/dt_coloc.txt")
dt_multicis <- fread( "Data/Modified/dt_multicis.txt")
df_VEP <- fread( "Data/Modified/df_VEP.txt")
res_multicis <- fread( "Data/Modified/res_multicis.txt")


res_unicis <- merge(dt_wald, dt_unicis, by = c("exposure", "outcome", "id.exposure", "id.outcome"), all.x = TRUE)
res_unicis <- merge(res_unicis, dt_gene_region[,.(hgnc, UniProt, id, study, gene_region)], by.x = "id.exposure", by.y ="id")
res_unicis <- merge(res_unicis, df_index[,.(id,clean_variable_name, note)], by.x = "id.outcome", by.y = "id")
setnames(res_unicis, "note", "outcome_category")

k<- df_VEP[, any(is_altering_variant), by = "SNP"]
setnames(k, "V1", "is_altering_variant")
res_unicis <- merge(res_unicis, k, by.x = c("nsnp"), by.y = c("SNP"), all.x = TRUE)

#include the meta-analysis
calculate_meta_results <- function(k) {
  arguments <- k[,.(hgnc, id.outcome)] %>% distinct
  dat_meta <- map(1:arguments[,.N], function(i) {
    dat<-k[hgnc==arguments[i, hgnc] & id.outcome == arguments[i,id.outcome],]
    meta_model <- rma.uni(yi = b, sei=se, data = dat, method = "REML")
    res<- data.table(hgnc=arguments[i, hgnc], id.outcome = arguments[i,id.outcome],
                     b = meta_model$beta[1,], se = meta_model$se, pval = meta_model$pval,
                     cochranQ = meta_model$QE, cochranQ_p=meta_model$QEp)
    return(res)
  }) %>% rbindlist(., fill = TRUE)
  if("posprob_coloc_PPH4" %in% colnames(k)){
    colnom <- paste0("posprob_coloc_PPH", 0:4)
    k<- k[,  lapply(.SD, mean), by = c("hgnc","outcome", "clean_variable_name", "id.outcome"), .SDcols = colnom]
    dat_meta <- merge(dat_meta, distinct(k[,.(hgnc,outcome, clean_variable_name,id.outcome, posprob_coloc_PPH0, posprob_coloc_PPH1, posprob_coloc_PPH2, posprob_coloc_PPH3, posprob_coloc_PPH4)]), by = c("hgnc", "id.outcome")) 
  } else {
    dat_meta <- merge(dat_meta, distinct(k[,.(hgnc,outcome, clean_variable_name,id.outcome)]), by = c("hgnc", "id.outcome")) 
  }
  dat_meta[,exposure := ""]
  dat_meta[,study:="meta-analysis"]
  dat_meta[,lci:=b-1.96*se]
  dat_meta[,uci:=b+1.96*se]
  return(dat_meta)
}
dat_meta <- calculate_meta_results(res_unicis[study %in% c("ARIC", "deCODE", "FENLAND"), ])
res_unicis<- rbind(res_unicis, dat_meta, fill = TRUE)

fwrite(res_unicis, "Data/Modified/res_unicis.txt")

res_multicis <- merge(res_multicis, dt_gene_region[,.(hgnc, UniProt, id, study, gene_region)], by.x = "id.exposure", by.y ="id")
res_multicis <- merge(res_multicis, df_index[,.(id,clean_variable_name, note)], by.x = "id.outcome", by.y = "id")
setnames(res_multicis, c("beta.multi_cis", "se.multi_cis", "pval.multi_cis"), c("b", "se", "pval"))
dat_meta <- calculate_meta_results(res_multicis[study %in% c("ARIC", "deCODE", "FENLAND"), ])
res_multicis<- rbind(res_multicis, dat_meta, fill = TRUE)
fwrite(res_multicis, "Data/Modified/res_multicis_meta.txt")

####
dt_multicis[method %in% c("Inverse variance weighted","Weighted median", "Contamination mixture" ), causal_multicis:= all(pval<0.05) & (all(b<0)|all(b>0)) , by = c("id.exposure", "id.outcome")]
dt_multicis <- merge(dt_multicis, dt_gene_region[,.(id, hgnc)], by.x = "id.exposure", by.y = "id")

#####
# ####You want a protein with potential of curing a disease (causal) with as little side effect as possible ####
# k <- res_unicis[study %in% c("deCODE", "ARIC", "FENLAND"),all(posprob_coloc_PPH4>0.5) & all(pval < 0.05) &(all(b<0)|all(b>0)), by = c("hgnc","outcome")]
# setnames(k, "V1", "causal")
# res_unicis <- merge(res_unicis, k, by = c("hgnc", "outcome"))
# therapeutic_target <- res_unicis[causal==TRUE & is_altering_variant != TRUE, mean(b),c("hgnc", "outcome")]
# # hgnc_toselect<- therapeutic_target[, all(V1<0)|all(V1>0), by="hgnc"][V1==TRUE,]$hgnc
# # therapeutic_target <- therapeutic_target[hgnc %in% hgnc_toselect,]
# therapeutic_target[,mean_causaleffect := mean(V1) , by = "hgnc"]
# side_effect <- merge(res_unicis, therapeutic_target[,.(hgnc, mean_causaleffect)] %>% distinct, by = "hgnc")
# # side_effect[pval<0.05 & ((b<0)!=(mean_causaleffect<0)), mean(b), by = "hgnc"]
# 
# #Bonus if colocalise with IUCPQ Biobank
# k1 <- merge(res_unicis[posprob_coloc_PPH4 > 0.5 & study %in% "IUCPQ Biobank", .(hgnc, outcome, id.outcome)], therapeutic_target, by = c("hgnc","outcome"))
# 
# 
# #Bonus if multi sign prioritise
# k2 <- merge(dt_multicis[causal_multicis == TRUE, .(hgnc, outcome, id.outcome)], therapeutic_target, by = c("hgnc","outcome")) %>% distinct
# 
# k<-rbind(k1, k2)
# res_sign <- k[,.(hgnc, id.outcome, outcome)] %>% distinct
# fwrite(res_sign, "Data/Modified/res_sign.txt")

message("This script finished without errors")
