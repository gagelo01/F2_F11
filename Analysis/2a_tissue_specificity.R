#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(tictoc)

wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/F2_finngen/"
setwd(wd)


coagulation_gene <- fread("Data/Modified/bloodcoagulationpathway_genes.txt") #genes from https://maayanlab.cloud/Harmonizome/gene_set/Blood+coagulation/PANTHER+Pathways

dt_gene_region <- fread("Data/Modified/dt_gene_region.txt")

gene_list <- dt_gene_region[hgnc %in% coagulation_gene$Symbol, hgnc %>% unique]

dt_united <- GagnonMR::get_tpm_for_genes_on_all_tissues(gene_list = gene_list)

GagnonMR::heatmap_tissular_specificity(dt_united = dt_united)

ggsave(paste0("Results/heatmap_tissular_specificiy_coagulation" ,".png"),
       width = 714/72,
       height = 630/72,
       units = "in",
       scale = 1,
       device="png")

fwrite("Data/Modified/dt_tissue_specificity.txt")
