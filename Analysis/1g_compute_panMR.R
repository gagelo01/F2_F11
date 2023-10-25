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
inst_pan <- fread( "Data/Modified/inst_pan.txt")
all_outcome <- fread( "Data/Modified/all_outcome.txt")


#### pan analysis
arguments <- crossing(inst_pan$id.exposure, all_outcome$id.outcome)
colnames(arguments)<-c("first", "second")
harm_univariate <- pmap(arguments, function(first, second) {
  TwoSampleMR::harmonise_data(exposure_dat = inst_pan[id.exposure == first, ],
                              outcome_dat = all_outcome[id.outcome == second, ],
                              action = 1)}) %>% rbindlist(.,fill = TRUE)

arguments <- harm_univariate[,.(id.exposure, id.outcome)] %>% distinct(.) #This line is to prevent id.exposure without SNP
colnames(arguments)<-c("first", "second")
setDT(harm_univariate)
harm_univariate[, exposure_outcome :=paste0(exposure, "_", outcome)]
egger_intercept <- TwoSampleMR::mr_pleiotropy_test(harm_univariate)
setDT(harm_univariate)

harm_univariate <- TwoSampleMR::steiger_filtering(harm_univariate) %>% as.data.table(.)
harm_univariate <- harm_univariate[mr_keep==TRUE,]
fwrite(harm_univariate, "Data/Modified/harm_univariate.txt")
harm_univariate <- harm_univariate[steiger_dir == TRUE,]


dt_pan <- future_pmap(arguments, function(first, second) {
  GagnonMR::all_mr_methods(harm_univariate[id.exposure == first & id.outcome == second,],
                           short = FALSE,
                           skip_presso = FALSE)}, 
  .options = furrr_options(seed = TRUE)) %>%
  rbindlist(., fill = TRUE)

FandQ <- future_pmap(arguments, function(first, second) {
  x<- harm_univariate[id.exposure == first & id.outcome == second,]
  res <- TwoSampleMR::mr_heterogeneity(x)
  if(nrow(res)==0){res<-data.frame(exposure = x$exposure, outcome = x$outcome)}
  x <- TwoSampleMR::add_rsq(x)
  res$fstat<-GagnonMR::fstat_fromdat(list(x))
  res$rsq <- sum(x$rsq.exposure)
  setDT(x)
  res$N <- x[,.N]
  return(res)
}, .options = furrr_options(seed = TRUE)) %>% 
  rbindlist(.,fill = TRUE)

setDT(egger_intercept)

list_res<-list(FandQ= FandQ, dt_pan=dt_pan, egger_intercept=egger_intercept)

fwrite(list_res$FandQ, "Data/Modified/FandQ_univariate.txt")
fwrite(list_res$dt_pan, "Data/Modified/dt_pan.txt")
fwrite(list_res$egger_intercept, "Data/Modified/egger_intercept.txt")

message("This script finished without errors")
