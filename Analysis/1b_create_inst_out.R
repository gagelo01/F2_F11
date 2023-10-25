#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
library(tictoc)
library(furrr)

wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/F2_finngen"
setwd(wd)
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref<-"/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
x <- paste0("awk -F ' ' '{print $2}' ","/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs", ".bim")
dt_ref_snp <- fread(cmd = x, header = FALSE) #Those are SNP in 1000G with MAF>0.01
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")


options(future.globals.maxSize= 5e9)
plan(multicore, workers = 20, gc = TRUE)
dt_gene_region <- fread("/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/PWMR_anything/Data/Modified/dt_gene_region.txt")
dt_gene_region[study=="Atherosclerosis Risk in Communities", study := "ARIC"]
newrow <- fread( "Data/Modified/newrow_eqtlgtex.txt")
coagulation_gene <- fread("Data/Modified/bloodcoagulationpathway_genes.txt") #genes from https://maayanlab.cloud/Harmonizome/gene_set/Blood+coagulation/PANTHER+Pathways
#####Modify dt_gene_region to include liver eqtl  #####
k <- newrow[, .(id, trait, consortium)]
k<-separate(k, col = "trait", into = c("tissue", "hgnc"), remove = FALSE)
k<- merge(k, distinct(dt_gene_region[,.(hgnc, UniProt, gene_region)]), by = "hgnc")
setnames(k, "consortium","study")
k[,tissue := NULL]
dt_gene_region <- rbind(dt_gene_region, k)

k<- df_index[grepl("eqtl-4-", id), .(id, clean_variable_name, consortium, trait)]
setnames(k, c("clean_variable_name", "consortium"), c("hgnc","study"))
k<- merge(k, distinct(dt_gene_region[,.(hgnc, UniProt, gene_region)]), by = "hgnc")
dt_gene_region <- rbind(dt_gene_region, k)
dt_gene_region[grepl("eqtl-74", id), vcffile := paste0(wd, "/Data/Modified/Gtex_vcf/", id, "/", id, ".vcf.gz")]
dt_gene_region[!grepl("eqtl-74", id), vcffile := paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", id, "/", id, ".vcf.gz")]

fwrite(dt_gene_region, "Data/Modified/dt_gene_region.txt")

######change parameters#########
######
parameters <- GagnonMR::default_param()
parameters$path <- c(GagnonMR::default_param()$path, paste0(wd, "/Data/Modified/Gtex_vcf/"))
parameters$uni_cis_minpval <- 1
parameters$ldref <- ldref
parameters$snp_bim <- dt_ref_snp$V1 #I only include SNP in 1000G with maf > 0.01
parameters$multicis_clumping["clump_r2"]<-0.4
############Choose the gene you wish to include and the outcome you wish to include  ###########
gene_toinclude <-  dt_gene_region[hgnc %in% coagulation_gene$Symbol, hgnc %>% unique]
ID_mrbase_out <- c(NULL)
out_mrbase <- NULL#paste0("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/", ID_mrbase_out, "/", ID_mrbase_out, ".vcf.gz")
ID_server_out <- c(
  "dis-15-157", "dis-15-1516", "dis-15-2495", "dis-15-2844", "dis-15-1739",#various bleeding outcome with ncase < 1000 in Finngen
                   "dis-1-1", "trait-7-1",  #venout_thromboembolism, parental lifespant (years),
                     df_index[grepl("dis-14-", id) & population == "European",id]) #Ischemicstroke, Cardioembolic stroke
out_server <- paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ID_server_out, "/", ID_server_out, ".vcf.gz")

####create inst_pan################
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")

list_study_together <- list(c("deCODE", "FENLAND")) #ARIC and IUCPQ Biobank only has cis summary statistics
arguments<-crossing(index = 1:length(list_study_together), gene_toinclude)
list_vcffile_vec <- pmap(arguments, function(index, gene_toinclude) dt_gene_region[study %in% list_study_together[[index]] & hgnc %in% gene_toinclude, ]$vcffile)

inst_pan <- future_map(list_vcffile_vec, function(x)
  {GagnonMR:::intern_inst_pan(vcffile_exp = x, parameters = parameters)},
           .options = furrr_options(seed = TRUE)) %>% rbindlist(., fill = TRUE)

fwrite(inst_pan, "Data/Modified/inst_pan.txt")

#######create inst_unicis########
list_study_together <- list(c("deCODE", "FENLAND", "ARIC"), c("IUCPQ Biobank")) #ARIC only has cis summary statistics
# list_study_together <- list(c("deCODE", "FENLAND", "ARIC", "GTEX", "IUCPQ Biobank")) #ARIC only has cis summary statistics

arguments<-crossing(index = 1:length(list_study_together), gene_toinclude)
list_vcffile_vec <- pmap(arguments, function(index, gene_toinclude)
  dt_gene_region[study %in% list_study_together[[index]] & hgnc %in% gene_toinclude, .(vcffile, gene_region)])

unicis_inst_from_vcffile_vec <- function(vcffile_vec,
                         gene_region,
                         parameters = default_param()
) {

  dat_tsmr <- GagnonMR:::intern_vcfpath_to_TwoSampelMR_region(vcf = vcffile_vec,
                                                              chrompos = gene_region, parameters = parameters)
  dt_null <- GagnonMR:::intern_dt_null(vcffile_exp = vcffile_vec,
                                       vcffile_out = NULL, parameters = parameters)
  rsid <- dat_tsmr[dat_tsmr$pval.exposure < parameters$uni_cis_minpval, ]$SNP
  if (!is.null(parameters$snp_bim)) {
    rsid <- rsid[rsid %in% parameters$snp_bim]
  }
  dat_ld <- dat_tsmr[dat_tsmr$SNP %in% rsid, ]
  if (dim(dat_ld)[1] == 0) {
    return(dt_null)
  }
  dat_ld <- dat_ld[, .SD[which.min(pval.exposure)], by = "id.exposure"]
  return(dat_ld)
}

inst_unicis <- future_map(list_vcffile_vec, function(x)
{unicis_inst_from_vcffile_vec(vcffile_vec = x$vcffile, gene_region = unique(x$gene_region), parameters = parameters)},
.options = furrr_options(seed = TRUE)) %>% rbindlist(., fill = TRUE)
setDT(inst_unicis)

fwrite(inst_unicis, "Data/Modified/inst_unicis.txt")

#####Create inst_multicis#######
inst_multicis_from_vcffile_vec <- function (vcffile_vec, chrompos, parameters = default_param()) {
  message("initializing multicis mr")
  stopifnot(gwasvcf::check_bcftools() & gwasvcf::check_plink())
  dat_tsmr <- GagnonMR:::intern_vcfpath_to_TwoSampelMR_region(vcf = vcffile_vec, 
                                                              chrompos = chrompos, parameters = parameters)
  dt_null <- GagnonMR:::intern_dt_null(vcffile_exp = vcffile_vec, vcffile_out = NULL, 
                            parameters = parameters)
  rsid <- dat_tsmr[dat_tsmr$pval.exposure < parameters$multicis_clumping["clump_p"], ]$SNP
  if (!is.null(parameters$snp_bim)) {
    rsid <- rsid[rsid %in% parameters$snp_bim]
  }
  dat_tsmr <- dat_tsmr[dat_tsmr$SNP %in% rsid, ]
  if (dim(dat_tsmr)[1] == 0) {
    return(dt_null)
  }
  dat_tsmr <- GagnonMR::clump_data_local(dat_tsmr, ldref = parameters$ldref, 
                                         clump_r2 = parameters$multicis_clumping["clump_r2"], 
                                         clump_kb = parameters$multicis_clumping["clump_kb"], 
                                         clump_p = 1)
  return(dat_tsmr)
}

list_study_together <- list("deCODE","FENLAND", "ARIC", "IUCPQ Biobank")
# arguments<-crossing(index = 1:length(list_study_together), gene_toinclude)
# list_vcffile_vec <- pmap(arguments, function(index, gene_toinclude)
#   dt_gene_region[study %in% list_study_together[[index]] & hgnc %in% gene_toinclude, .(vcffile, gene_region)])
k <- dt_gene_region[study%in%list_study_together & hgnc %in% gene_toinclude, .(vcffile, gene_region) ] #otherwise to much instrument
list_vcffile_vec <- split(k, 1:k[,.N])
inst_multicis <- future_map(list_vcffile_vec, function(x)
{inst_multicis_from_vcffile_vec(vcffile_vec = x$vcffile, chrompos = unique(x$gene_region), parameters = parameters)},
.options = furrr_options(seed = TRUE)) %>% rbindlist(., fill = TRUE)
setDT(inst_multicis)

fwrite(inst_multicis, "Data/Modified/inst_multicis.txt")

#####create out#####
inst_pan <- fread( "Data/Modified/inst_pan.txt")
inst_unicis <- fread( "Data/Modified/inst_unicis.txt")
inst_multicis <- fread( "Data/Modified/inst_multicis.txt")

tic()
out <- future_map(as.list(c(out_mrbase, out_server)), function(x, rsiid = unique(c(inst_unicis[!is.na(SNP),]$SNP, inst_pan[!is.na(SNP),]$SNP, inst_multicis[!is.na(SNP),]$SNP))) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  outr<-GagnonMR::extract_outcome_variant(snps = rsiid, outcomes = x, rsq = 0.6, parameters = parameters)
  return(outr)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)
toc()

fwrite(out, "Data/Modified/all_outcome.txt" )

message("This script finished without errors")
