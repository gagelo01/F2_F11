#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
library(furrr)
library(tictoc)
wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/F2_finngen/"
setwd(wd)

#####load objects####
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
coagulation_gene <- fread("Data/Raw/bloodcoagulationpatway_genes.txt") #genes from https://maayanlab.cloud/Harmonizome/gene_set/Blood+coagulation/PANTHER+Pathways
coagulation_gene <- rbind(coagulation_gene, data.table(Symbol = "F11", Name = "coagulation factor XI"))
fwrite(coagulation_gene, "Data/Modified/bloodcoagulationpathway_genes.txt")
dt_gene_region <- fread("/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/PWMR_anything/Data/Modified/dt_gene_region.txt")
gencode <- fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt")
ref_data = "/mnt/sda/couchr02/Eloi/b37_b38/ref_75MSNPS_b37tob38_ALFA_EUR.txt.gz"

######The arguments you need to change#######
vec_gene <-  dt_gene_region[hgnc %in% coagulation_gene$Symbol, hgnc %>% unique]
vec_tissue <-"Liver"
mywindow = 2e6

#########get dt_eqtl#########
gwasvcf::set_bcftools()
gwasvcf::set_plink()
options(future.globals.maxSize= 5e9)
plan(multicore, workers = 20)
dt_eqtl <-  future_map(as.list(vec_gene), function(x) {
  GagnonMR::get_eQTL(tissue = vec_tissue, gene = x, mywindow = mywindow)  }, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)


fwrite(dt_eqtl, "Data/Modified/dt_eqtl.txt")

##############make gtex as vcf###################
info <- dt_eqtl[, .(exposure, sample_size, gene_name)] %>% distinct(.)
dt_chrompos <- dt_eqtl[, paste0(unique(chr), ":", min(pos), "-", max(pos)), by = "exposure"]
setnames(dt_chrompos, "V1", "chrompos")
info <- merge(info, dt_chrompos, by = "exposure")
id_prefix <- "eqtl-74-"

newrow <- data.table(id = paste0(id_prefix, 1:info[,.N]), trait = info$exposure, group_name = "public",
                     year = 2020, author = "The GTEX Consortium",consortium = "GTEX",
                     sex = "Males and Females", population = "European",
                     initial_build = "HG38/GRCh38", unit ="SD", nsnp = "dummy", sample_size = info$sample_size,
                     category = "eqtl", pmid = 32913098, ncase = NA,sd = 1, note = info$chrompos, ncontrol = NA,
                     adjustments = "Adjusted for the two first PEER factors (probabilistic estimation of expression
                     residuals), the two first genetic principal components of ancestry, the platform used and sex",
                     clean_variable_name = info$exposure)


df_index <- df_index[!grepl(id_prefix, id), ]
df_index <- rbind(df_index, newrow)
# Delete Directory
map(newrow$id, function(x) {
  dir <- paste0(wd, "/Data/Modified/Gtex_vcf/", x)
  if (file.exists(dir)) {
    unlink(dir,recursive = TRUE)
    message(paste0("Directory ", x," has been deleted"))
  }
})
#create vcf
dir.create("Data/Modified/Gtex_vcf")


create_my_vcf <-function(i, newrow, df_index) {
  message("Creating VCF for ", paste0(newrow[i, id]))

  y <- paste0("tabix -h ",ref_data," ", newrow[i, ]$note)
  myref <- fread(cmd=y)
  setnames(myref, c("#chr", "EUR_af", "ref", "alt"), c("chr", "eur_eaf", "other_all", "effect_all"))
  myref[, eur_eaf := eur_eaf %>% ifelse(.==0,0.001,. ) %>% ifelse(.==1, 0.999, .)]

  dt_eqtl <- dt_eqtl[!(is.na(se)|is.na(eaf)),]

  GagnonMR::formattovcf_createindex2(all_out = dt_eqtl,
                           snp_col = NULL,
                           outcome_name = newrow[i, trait],
                           beta_col = "beta",
                           se_col = "se",
                           pval_col = "p_val",
                           eaf_col = "eaf",
                           maf_col = NULL,
                           effect_allele_col = "effect_allele",
                           other_allele_col =  "other_allele",
                           ncase_col = NULL,
                           ncontrol_col = NULL,
                           samplesize_col = "sample_size",
                           chr_col = "chr",
                           pos_b37_col = NULL,
                           pos_b38_col = "pos",
                           units = "SD",
                           traduction = myref,
                           out_wd = paste0(wd, "/Data/Modified/Gtex_vcf"),
                           df_index = df_index,
                           group_name = newrow[i, group_name],
                           year = newrow[i, year],
                           author = newrow[i, author],
                           consortium = newrow[i, consortium],
                           sex = newrow[i, sex],
                           population = newrow[i, population],
                           initial_build = newrow[i, initial_build],
                           category = newrow[i, category],
                           pmid = newrow[i, pmid],
                           note = newrow[i, note],
                           ID = newrow[i,id])
}


tictoc::tic()
future_map(as.list(c(1:newrow[,.N])), function(x) {
  create_my_vcf(i = x , newrow = newrow, df_index = df_index)}, .options = furrr_options(seed = TRUE))
tictoc::toc()
fwrite(newrow, "Data/Modified/newrow_eqtlgtex.txt")

message("This script finished without errors")
