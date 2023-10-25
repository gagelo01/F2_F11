#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)
setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/F2_finngen/")
dt_wald <- fread( "Data/Modified/dt_wald.txt")
dt_wald <- dt_wald[,.(id.exposure, id.outcome, outcome, exposure, method, nsnp, b, se,
                                    pval,type_of_test, lci, uci, fstat, rsq, steiger_dir, steiger_pval)]
dt_gene_region <- fread("Data/Modified/dt_gene_region.txt")
inst_unicis <- fread( "Data/Modified/inst_unicis.txt")
ref_data <- "/mnt/sda/couchr02/Eloi/b37_b38/ref_78MSNPS_b37tob38_ALFA_EUR.txt.gz"

gwasvcf::set_bcftools()
gwasvcf::set_plink()

mytemp <- paste0(tempfile(), ".txt")
k <- inst_unicis[SNP %in% dt_wald$nsnp, .(chr.exposure, pos.exposure)] %>% distinct(.)
k[,pos:=pos.exposure]
fwrite(k, mytemp, col.names = FALSE, sep = "\t")

x<- paste0("tabix -h ", ref_data, " -R ", mytemp)
k <- fread(cmd = x)
all_out_vcf <- gwasvcf::create_vcf(chrom= k[,`#chr`], pos= k[, pos_b37],	
                                   snp= k[, rsid], nea= k[,ref_b37],	
                                   ea= k[,alt_b37], name="1000g")
VariantAnnotation::writeVcf(all_out_vcf, file = "Data/Modified/VEPinput.vcf")

########

system2(command = "sh", args = "Analysis/1f_vep.sh")

#
# vepoutput <-fread("Data/Modified/VEPoutput.txt", skip = 39)
# test <- all_out_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(.) %>% as.data.table
vepoutput <-fread("Data/Modified/VEPoutput.txt", skip = 25)
# test <- all_out_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(.) %>% as.data.table
setnames(vepoutput, "#Uploaded_variation", "Uploaded_variation")
vepoutput <- separate(vepoutput, "Consequence", sep = ",", into = paste0("consequence", 1:4))
vepoutput <- vepoutput[, unique(c(consequence1, consequence2, consequence3, consequence4)), by =c("Uploaded_variation", "SYMBOL")]
setnames(vepoutput, "V1", "consequence")
vepoutput<-vepoutput[!is.na(consequence)]
vepoutput[, consequences := paste(consequence, collapse = ","), by = c("Uploaded_variation", "SYMBOL")]
setnames(vepoutput, "Uploaded_variation", "SNP")
df_VEP <- distinct(vepoutput[, c("SNP", "SYMBOL", "consequences")])
toMatch<- c("missense_variant", "stop_gained", "stop_lost", "start_gained", "start_lost", "frameshift")
df_VEP[ , is_altering_variant := grepl(paste(toMatch,collapse="|"), consequences)]

fwrite(df_VEP, "Data/Modified/df_VEP.txt")

message("This script finished without errors")

