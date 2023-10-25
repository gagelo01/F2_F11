#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(metafor)

wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/F2_finngen"
setwd(wd)
source("Analysis/my_phenotypePlot_function.R")
res_unicis <- fread ( "Data/Modified/res_unicis.txt")
dt_wald <- fread("Data/Modified/dt_wald.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
res_multicis <- fread("Data/Modified/res_multicis_meta.txt")
dt_pan <-fread("Data/Modified/dt_pan.txt")
res_hypr <- fread("Data/Modified/res_hypr.txt")
inst_unicis <- fread( "Data/Modified/inst_unicis.txt")
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
######
return_format_data<-function(data) {
k <- data[, paste0("OR = ", format(round(exp(b), digits = 2), nsmall = 2), ", 95% CI=", format(round(exp(lci), digits = 2), nsmall = 2), "-",  format(round(exp(uci), digits = 2), nsmall = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))]
names(k) <- data[,paste0(exposure," on ", outcome)]
return(k)
  }
return_format_data_noexp <-function(data) {
  return(data[, paste0(format(round(b, digits = 2), nsmall = 2) , " 95% CI=", format(round(lci, digits = 2), nsmall = 2), "-",  format(round(uci, digits = 2), nsmall = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))])
}
return_format_fstat <-function(data) {
  return(data[, paste0(N, " SNPs (r2 = ", round(rsq*100, digits =2), "%; F-statistics = ",  round(fstat, digits = 0), ")")])
}

####Abstract####
res_unicis$hgnc %>% unique %>% length()
df_index[id %in% c(res_unicis$id.exposure, res_unicis$id.outcome), .(pmid, consortium, ncase, ncontrol, sample_size)] %>% distinct
dat_meta <- res_unicis[ hgnc%in%c("F11", "F2")& study == "deCODE" & exposure%in%c("2190_55", "5316_54"),]
dat_meta[,b:=b*-1]
dat_meta[,lci:=b-1.96*se]
dat_meta[,uci:=b+1.96*se]

dat_meta[id.outcome%in%c("dis-1-1", "dis-14-8", "dis-15-157") & hgnc%in%c("F2", "F11"),] %>% return_format_data(.)

#####Intro######
res_unicis$hgnc %>% unique %>% length()


####Results#####
#Vte GWAS
x <- gwasvcf::query_gwas("/mnt/sda/gagelo01/Vcffile/Server_vcf/dis-1-1/dis-1-1.vcf.gz", pval = 5e-8)
x<- x%>%gwasvcf::vcf_to_tibble(.,)
setDT(x)
x[,length(unique(rsid))]
k<- get_inst(vcffile = "/mnt/sda/gagelo01/Vcffile/Server_vcf/dis-1-1/dis-1-1.vcf.gz", r2 = 0.01, kb = 1000)
k[,.N]
#uni-cis MR
res_unicis[study %in% c("deCODE", "FENLAND", "ARIC"), max(fstat), by = "hgnc"][V1<10,]
res_unicis[study %in% c("deCODE", "FENLAND", "ARIC"), mean(fstat)]
res_unicis[study == "deCODE" & is_altering_variant==TRUE,unique(hgnc)]
res_unicis[study == "deCODE", unique(hgnc)%>%length] - res_unicis[study == "deCODE" & is_altering_variant==TRUE,length(unique(hgnc))]
res_unicis[steiger_dir==FALSE,]

res_unicis[study=="meta-analysis" & id.outcome %in% c("dis-1-1", "dis-14-7"), ][order(-abs(b)), .(hgnc, outcome, b, posprob_coloc_PPH4)]

dat_meta[id.outcome=="dis-1-1"&hgnc=="F2"] %>% return_format_data(.)
dat_meta[id.outcome=="dis-14-8"&hgnc=="F2"] %>% return_format_data(.)
dat_meta[id.outcome=="dis-15-157"&hgnc=="F2"] %>% return_format_data(.)
dat_meta[id.outcome=="dis-1-1"&hgnc=="F11"] %>% return_format_data(.)
dat_meta[id.outcome=="dis-14-8"&hgnc=="F11"] %>% return_format_data(.)
dat_meta[id.outcome=="dis-15-157"&hgnc=="F11"] %>% return_format_data(.)

dat_meta2<-res_unicis[hgnc=="SERPINC1" & study == "deCODE",] 
dat_meta2[id.outcome=="dis-1-1"] %>% return_format_data(.)
dat_meta2[id.outcome=="dis-14-7"] %>% return_format_data(.)
dat_meta2[id.outcome=="dis-15-157"] %>% return_format_data(.)
res_unicis[id.outcome%in%c("dis-14-7", "dis-14-7", "dis-1-1", "dis-15-157")&hgnc%in%c("F11", "F2")&study=="meta-analysis", .(hgnc, id.outcome, posprob_coloc_PPH4)] 

res_unicis[id.outcome%in%c("dis-14-7", "dis-1-1", "dis-15-157")&hgnc%in%c("F11", "F2")&study=="meta-analysis", min(cochranQ_p)] 

#Multi-cis MR
res_multicis[id.outcome%in%c("dis-14-7", "dis-1-1", "dis-15-157")&hgnc%in%c("F11", "F2")&study=="meta-analysis", min(cochranQ_p)] 

meta_multicis<-res_multicis[hgnc %in% c("F2", "F11") & id.outcome%in%c("dis-14-7", "dis-1-1", "dis-15-157")
                            & study == "deCODE" & exposure%in%c("2190_55", "5316_54"), ]
# meta_multicis[, min(cochranQ_p)]
meta_multicis[,b:=b*-1]
meta_multicis[,lci:=b-1.96*se]
meta_multicis[,uci:=b+1.96*se]
meta_multicis%>% return_format_data(.)

#Pan MR
dt_pan<-merge(dt_pan, dt_gene_region, by.x = "id.exposure", by.y = "id")
res_pan <- dt_pan[hgnc%in%c("F11", "F2")&
         id.outcome%in%c("dis-14-7", "dis-1-1", "dis-15-157") & 
           study %in% c("deCODE", "FENLAND") & exposure%in%c("2190_55", "5316_54", "fenland2190_55", "fenland5316_54"),]
res_pan[,b:=b*-1]
res_pan[,lci:=b-1.96*se]
res_pan[,uci:=b+1.96*se]
res_pan[study == "deCODE" & method=="Inverse variance weighted"]  %>% return_format_data(.)
themet<-c("Inverse variance weighted", "Weighted median", "MR-PRESSO (outlier-corrected)", 
          "Contamination mixture", "MR Egger", "Weighted mode")
res_pan[method %in% themet, (all(b<0)|all(b>0))&all(pval<0.05), by=c("id.exposure", "id.outcome", "hgnc", "exposure") ]
dt_pan[hgnc%in%c("F11", "F2")&
         id.outcome%in%c("dis-14-7", "dis-1-1", "dis-15-157") & 
         study == "deCODE" & method %in% themet, (all(b<0)|all(b>0))&all(pval<0.05), by=c("hgnc","exposure", "id.outcome")]

dt_pan[exposure=="5316_54"&id.outcome %in%c("dis-14-7", "dis-1-1") & pval>0.05,]
#hepatic gene expression

k<- res_unicis[study %in% "IUCPQ Biobank" & hgnc %in% c("F2", "F11")& id.outcome%in%c("dis-14-7", "dis-1-1","dis-15-2844", "dis-15-157", "trait-7-1"), ]
inst_unicis[id.exposure %in% k$id.exposure, pval.exposure %>% formatC(., format = "e", digits = 1)]
k[,(min(rsq))%>%round(.,digits = 2) ,by = "hgnc"]
k[,min(fstat)%>%round(.,digits = 0) ,by = "hgnc"]
k[,b:=b*-1]
k[,lci:=b-1.96*se]
k[,uci:=b+1.96*se]
k %>% return_format_data()
k[id.outcome == "trait-7-1", ]%>% return_format_data_noexp()
k[, .(hgnc, id.outcome, posprob_coloc_PPH3,posprob_coloc_PPH4)]
res_hypr


#hyprcoloc
studytoinclude<-c("FENLAND", "deCODE", "ARIC")
res_coloc[study_exposure %in% studytoinclude & study_outcome %in% studytoinclude, mean(posprob_coloc_PPH4), by = "hgnc"]

res_coloc[(study_exposure %in% studytoinclude & study_outcome %in% "IUCPQ Biobank") |
            (study_exposure %in% "IUCPQ Biobank" & study_outcome %in% studytoinclude) , mean(posprob_coloc_PPH4), by = "hgnc"]

res_coloc[(study_exposure %in% studytoinclude & study_outcome %in% "IUCPQ Biobank") |
            (study_exposure %in% "IUCPQ Biobank" & study_outcome %in% studytoinclude) , mean(posprob_coloc_PPH3), by = "hgnc"] #So very clearly F11 protein and gene expression do not share the same causal SNP


#####Methods######
res_unicis[, length(unique(hgnc)), by = "study"]
res_unicis[, length(unique(exposure)), by = "study"]


df_index[grepl("dis-14-", id) & id %in% dt_pan$id.outcome, .(trait,ncase, ncontrol)]
df_index[grepl("dis-15-", id) & id %in% dt_pan$id.outcome, .(trait,clean_variable_name)]
dt_trad <- fread("/mnt/sda/gobemi01/PheWAS/PheWAS_MR_Lpa/data/FinnGen_r7_infofile.csv")
dt_trad[phenocode=="BLEEDING",]
library(readxl)
info_finngen <- readxl::read_xlsx("Data/Raw/FINNGEN_ENDPOINTS_DF7_Final_2021-03-05_public.xlsx")
setDT(info_finngen)
info_finngen[NAME == "BLEEDING",]$INCLUDE
info_finngen[NAME %in% df_index[grepl("dis-15-", id) & id %in% dt_pan$id.outcome, ]$trait,.(NAME, HD_ICD_10, HD_ICD_10_EXCL)]
info_finngen[NAME == "I9_INTRACRA",]$INCLUDE

