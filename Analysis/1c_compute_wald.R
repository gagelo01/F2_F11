#!/usr/bin/env Rscript
library(TwoSampleMR)
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
inst_unicis <- fread( "Data/Modified/inst_unicis.txt")
all_outcome <- fread( "Data/Modified/all_outcome.txt")

#### dt_uni_cis ####

dt_wald <- future_map(as.list(inst_unicis$id.exposure %>% unique),
                             function(x) {
                               wald_quick(inst_all_sign_clump_arg = inst_unicis,
                                          all_outcome_arg = all_outcome, id_exposure_name = x)
                             }, .options = furrr_options(seed = TRUE)) %>% rbindlist(., fill = TRUE)

fwrite(dt_wald, "Data/Modified/dt_wald.txt")

message("This script finished without errors")
