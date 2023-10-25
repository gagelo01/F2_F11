#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(furrr)
wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/F2_finngen"
setwd(wd)
source("Analysis/func_F2.R")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
options(future.globals.maxSize= 5e9)
plan(multicore, workers = 20, gc = TRUE)
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
inst_multicis <- fread( "Data/Modified/inst_multicis.txt")
all_outcome <- fread( "Data/Modified/all_outcome.txt")

####perform multi-cis MR####
parameters <- GagnonMR::default_param()
k <- inst_multicis[!is.na(beta.exposure), SNP, by = "id.exposure"] %>%
  split(., .$id.exposure) 
list_snp <- map(k, function(x) x$SNP)
list_ldmat <- future_map(list_snp, function(x) {
ldmat <- ieugwasr::ld_matrix_local(unique(x), 
                                   plink_bin = genetics.binaRies::get_plink_binary(), bfile = parameters$ldref)
return(ldmat)
})

arguments <- crossing(inst_multicis$id.exposure, all_outcome$id.outcome)
colnames(arguments)<-c("first", "second")

multicis_get <- function(first, second) {
harm <-  TwoSampleMR::harmonise_data(exposure_dat = inst_multicis[id.exposure == first, ],
                                                     outcome_dat = all_outcome[id.outcome == second, ],
                                                     action = 1)

if(nrow(harm)<2){return(NULL)}
ldmat <- list_ldmat[[first]]
x<-harm %>% as.data.table(.)
mriobj <- MendelianRandomization::mr_input(bx = x$beta.exposure, 
                                           bxse = x$se.exposure, by = x$beta.outcome, byse = x$se.outcome, 
                                           exposure = x$exposure[1], outcome = x$outcome[1], 
                                           snps = x$SNP, effect_allele = x$effect_allele.exposure, 
                                           other_allele = x$other_allele.exposure, eaf = x$eaf.exposure)
k <- x[, c(paste(SNP, effect_allele.exposure, other_allele.exposure, 
                 sep = "_"), paste(SNP, other_allele.exposure, 
                                   effect_allele.exposure, sep = "_"))]
mriobj@correlation <- ldmat[k[k %in% rownames(ldmat)], 
                            k[k %in% colnames(ldmat)]]
IVWobject <- MendelianRandomization::mr_ivw(mriobj, 
                                            model = "default", correl = TRUE, distribution = "normal", 
                                            alpha = 0.05)
res_multicis <- cbind(x[1, c("id.exposure", "id.outcome", 
                             "exposure", "outcome")], data.frame(beta.multi_cis = IVWobject@Estimate, 
                                                                 se.multi_cis = IVWobject@StdError, pval.multi_cis = IVWobject@Pvalue, 
                                                                 cochranQ.multi_cis = IVWobject@Heter.Stat[1], 
                                                                 cochranQpval.multi_cis = IVWobject@Heter.Stat[2], 
                                                                 nsnp.multi_cis = nrow(x)))
return(res_multicis)
}

multicis_get_safely <- safely(multicis_get)
res_multicis <- pmap(arguments, function(first, second) {
  multicis_get_safely(first=first, second=second)
}) 

res_multicis <- lapply(res_multicis, function(x) x$result)%>% rbindlist(.,fill=TRUE)

fwrite(res_multicis, "Data/Modified/res_multicis.txt")

message("This script finished without errors")
