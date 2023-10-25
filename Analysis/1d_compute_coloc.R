#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(furrr)

wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/F2_finngen"
setwd(wd)
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
options(future.globals.maxSize= 5e9)
plan(multicore, workers = 20, gc = TRUE)
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
dt_wald <- fread( "Data/Modified/dt_wald.txt")
parameters <- GagnonMR::default_param()
parameters$path <- c(GagnonMR::default_param()$path, paste0(wd, "/Data/Modified/Gtex_vcf/"))

#####
arguments <- merge(dt_wald[, .(id.exposure, id.outcome)], dt_gene_region[,.(id, vcffile, gene_region)], by.x = "id.exposure", by.y = "id")
# arguments <- merge(dt_wald[pval<0.05, .(id.exposure, id.outcome)], dt_gene_region[,.(id, vcffile, gene_region)], by.x = "id.exposure", by.y = "id")
arguments[,vcffile_out := paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", id.outcome, "/", id.outcome, ".vcf.gz")]
setnames( arguments, c("vcffile", "gene_region"), c("vcffile_exp", "chrompos"))
arguments <- arguments[,.(vcffile_exp, vcffile_out, chrompos)]

dt_coloc <- furrr::future_pmap(arguments, function(vcffile_exp, vcffile_out, chrompos)
                               GagnonMR::get_coloc(vcffile_exp= vcffile_exp, 
                                                   vcffile_out = vcffile_out,
                                                   chrompos = chrompos,
                                                   parameters = parameters)
  ,.options = furrr_options(seed = TRUE)) %>% rbindlist(., fill = TRUE)

fwrite(dt_coloc, "Data/Modified/dt_coloc.txt")
message("This script finished without errors")
