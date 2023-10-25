#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(furrr)
library(MetBrewer)
library(gassocplot)
library(metafor)
library(gassocplot)
library(forestploter)

wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/F2_finngen"
setwd(wd)
source("Analysis/my_phenotypePlot_function.R")
res_unicis <- fread ( "Data/Modified/res_unicis.txt")
res_sign <- fread( "Data/Modified/res_sign.txt")
res_unicis[clean_variable_name=="", clean_variable_name:=outcome]
res_multicis <- fread("Data/Modified/res_multicis_meta.txt")
dt_pan <-fread("Data/Modified/dt_pan.txt")
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
efficacy_outcome <- c("dis-1-1", "dis-14-7", "dis-14-8","dis-14-9", "dis-14-10")
safety_outcome <- c("dis-15-157", "dis-15-1516", "dis-15-1739", "dis-15-2844", "trait-7-1")

#balloon plot
#Fixed meta-analysis for each gene
from_resunicis_to_datballoon <- function(dat) {
  dat_balloon <- dat
  dat_balloon <- dat_balloon[,.(id.outcome, hgnc, clean_variable_name, b, pval, posprob_coloc_PPH4)]
  setnames(dat_balloon, c("hgnc", "clean_variable_name", "b", "pval"), 
           c("exposure", "outcome", "z_score", "pval"))
  dat_balloon <- distinct(dat_balloon)
  dat_balloon[, Category := ifelse(id.outcome %in%efficacy_outcome, "Efficacy outcomes", "Safety outcomes")]
  dat_balloon[,shape_point:=ifelse(pval<0.05,"circle", "cross")]
  dat_balloon[, Category := factor(Category, levels = c("Efficacy outcomes", "Safety outcomes"))]
  dat_balloon <- dat_balloon[order(Category, outcome, exposure),]
  dat_balloon[, outcome := factor(outcome, levels = unique(outcome)),]
  dat_balloon[,log10_pval:=posprob_coloc_PPH4]
  dat_balloon[, id.outcome := factor(id.outcome, levels = c(efficacy_outcome, safety_outcome))]
  dat_balloon[,outcome := factor(outcome, levels = dat_balloon[order(id.outcome), unique(outcome)])]
}


dat_balloon <- from_resunicis_to_datballoon(res_unicis[id.outcome%in%c(efficacy_outcome, safety_outcome) &
                                                         study == "meta-analysis"])
plot_balloon2(dat_balloon, subdivide_by_category = TRUE, name_pval = "PPH4", name_zcore = "Beta")

ggsave(filename = paste0("Results/balloonplot_bloodprotein", ".png"),
       width = 423/72,height = 642/72,units="in",scale=1, device = "png")

dat_balloon <- from_resunicis_to_datballoon(res_unicis[id.outcome%in%c(efficacy_outcome, safety_outcome) &
                                                         study%in%c("IUCPQ Biobank"),])
plot_balloon2(dat_balloon, subdivide_by_category = TRUE, name_pval = "PPH4", name_zcore = "Beta")

ggsave(filename = paste0("Results/balloonplot_hepaticgene", ".png"),
       width = 423/72,height = 642/72,units="in",scale=1, device = "png")

#####scatter plot pan#####
dt_pan <-fread("Data/Modified/dt_pan.txt")
dt_pan[method %in% c("Inverse variance weighted","Weighted median", "Contamination mixture", "MR-PRESSO (outlier-corrected)"), causal_multicis:= all(pval<0.05) & (all(b<0)|all(b>0)) , by = c("id.exposure", "id.outcome")]
harm_univariate <- fread( "Data/Modified/harm_univariate.txt")
cleanify <- function(dt) {
  dt<- merge(dt, df_index[, .(id,clean_variable_name)], by.x = "id.outcome", by.y = "id")
  dt <- merge(dt, dt_gene_region[,.(id, hgnc, study)], by.x = "id.exposure", by.y = "id")
  dt[, label := paste0(hgnc, "_", study)]
  return(dt)
}
list_pan <- map(list(dt_pan, harm_univariate), cleanify)
outcome_to_include <- c("Ischemic stroke", "Venous thromboembolism", "Bleeding")

k2 <-list_pan[[1]][hgnc%in%c("F11", "F2")&
                     clean_variable_name %in% outcome_to_include & 
                     study %in% c("deCODE", "FENLAND") & exposure%in%c("2190_55", "5316_54", "fenland2190_55", "fenland5316_54")]
k2 <- k2[method%in%c("MR-PRESSO (outlier-corrected)", "Contamination mixture",
                     "Inverse variance weighted", "Weighted median", 
                     "MR Egger", "Weighted mode")]
k2 <- k2[order(hgnc, study),]
GagnonMR::my_mr_scatter_plot(mr_results = k2,
                             dat = list_pan[[2]][exposure%in%unique(k2$exposure) & outcome %in% unique(k2$outcome), ][order(hgnc, study),],
                             equation_facet_grid = "label~clean_variable_name", legend.position = "top") +
  xlab("SNP effect on blood protein levels") +
  ylab("SNP effect on diseases") +
  guides(colour=guide_legend(nrow=2,byrow=TRUE))

ggsave(filename = paste0("Results/F2F11scatter_plot.png"),
       width = 700/72,height = 700/72,units="in",scale=1, device = "png")


# #### hyprcoloc plot #####
# res_hypr <- fread("Data/Modified/res_hypr.txt")
# exposure_hypr <- fread( "Data/Modified/exposure_hyp.txt")
# 
# dir.create("Results/Hypr_plot")
# exposure_vec <- c("F11", "F2")
# df_aligned <- exposure_hypr[gene.exposure %in% exposure_vec, ] %>% distinct
# studyorder<-c("IUCPQ Biobank"=1, "ARIC"=3, "deCODE"=4, "FENLAND"=5, "FinnGen"=6, "dis-1-1"=7, "dis-14-7"=8, "dis-14-8"=9)
# df_aligned[, study_order := studyorder[study]]
# 
# res_hypr1 <- res_hypr[gene_investigated %in% exposure_vec & iteration == 1, ]
# # argument <- df_aligned[,.(gene.exposure, id.exposure, study)] %>% distinct
# # argument <- argument[, .SD[1], by = c("study", "gene.exposure")]
# # argument <- split(argument,list(argument$gene.exposure))
# argument <- res_hypr1[,.(allexposures, gene_investigated)]
# 
# for(i in 1:length(argument)) {
#   k<- argument[i,]
#   theidexposure <- str_split(k$allexposures, ":")
#   k1<-df_aligned[gene.exposure %in% k$gene_investigated & id.exposure %in% theidexposure[[1]],]
#   k2<-res_hypr1[gene_investigated %in% unique(k$gene_investigated) & allexposures %in% k$allexposures,]
#   # if(k$gene.exposure[1]=="F2"){k2$candidate_snp<-"rs3136516"}
# 
#   k1[,exposure := paste0(study, "_", clean_variable_name)]
#   A <- stack_assoc_plot_wrapper(df_aligned = k1,
#                                 res_hypr1 = k2,
#                                 ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs",
#                                 traits_inorder = k1[order(study_order), unique(exposure)],
#                                 build = 37)
#   res <- sensitivity.plot_wrapper(df_aligned = k1,
#                                   traits_inorder = unique(k1$exposure))
#   B<-drawheatmap(res[[2]])
#   twopanel <-  cowplot::ggdraw() +
#     cowplot::draw_plot(ggplotify::as.ggplot(A) + theme(text = element_text(size = 0.4)), x = 0.08, y =0, width = .6, height = 1) +
#     cowplot::draw_plot(B, x = .66, y =0.1, width = .31, height = 0.7) +
#     cowplot::draw_plot_label(label = c("", ""), size = 25,
#                              x = c(0, 0.62), y = c(0.9, 0.9))
#   ggsave(plot = A, filename = paste0("Results/Hypr_plot/", "stackassoc_plot_",k[1]$gene.exposure, ".png"),
#          width = 790/72,height = 683/72,units="in",scale=1, device = "png")
#   saveRDS(object = twopanel, file = paste0("Results/Hypr_plot/", "twopanel_hypr_plot_",k[1]$gene.exposure, ".rds"))
#   ggsave(plot = twopanel, filename = paste0("Results/Hypr_plot/", "twopanel_hypr_plot_",k[1]$gene.exposure, ".png"),
#          width = 790/72,height = 683/72,units="in",scale=1, device = "png")
# }
# 
# 
# 
# message("This script finished without errors")

###Table for coloc
table_coloc <- res_unicis[hgnc %in% c("F2", "F11") & study %in% c("meta-analysis", "IUCPQ Biobank") & id.outcome %in% c(efficacy_outcome, safety_outcome),]
table_coloc[, Category := ifelse(id.outcome %in%efficacy_outcome, "Efficacy outcomes", "Safety outcomes")]
table_coloc[, exposure := gsub("SeqId_|fenland", "", exposure)]
table_coloc <- table_coloc[, .(Category, clean_variable_name, hgnc, study ,exposure, 
                               posprob_coloc_PPH0, posprob_coloc_PPH1, posprob_coloc_PPH2, posprob_coloc_PPH3,posprob_coloc_PPH4), ][order(Category, clean_variable_name, hgnc, study ,exposure)]
table_coloc[, omic := ifelse(study == "IUCPQ Biobank", "Hepatic gene levels", "Blood protein levels")]
my_theme<- ggplot2::theme(
  legend.title=element_blank(),
  panel.grid.major.y = element_line(size = 0.25, colour = "gray60"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.background = element_blank(),
  plot.margin = margin(t = 2, r = 0.5, b = 0.5, l = 0.5, "cm"),
  legend.position = "top",
  legend.text = element_text(
    color = "gray20",
    size = 10,
    margin = margin(l = 0.2, r = 0.2)
  ),
  legend.spacing.y = unit(0.1, 'cm'),
  legend.key = element_rect(fill = "transparent", colour = "transparent"),
  legend.key.size = unit(0.8, "cm"),
  axis.line = element_line(size = 0.5, colour = "gray20"),
  axis.ticks = element_line(size = 0.5, colour = "gray20"),
  axis.text.y = element_text(
    size = 10,
    colour = "gray20"
  ),
  axis.text.x = element_text(
    angle = 60,
    size = 8,
    hjust = 1,
    face = "plain",
    colour = "gray20"
  ))

ggplot(table_coloc, aes(x=clean_variable_name, y=posprob_coloc_PPH4, fill=omic))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  facet_grid(hgnc ~ Category, scales = "free_x") +
  # geom_text(aes(label=posprob_coloc_PPH4), vjust=1.6, color="white",
  #           position = position_dodge(0.9), size=3.5)+
  scale_fill_brewer(palette="Paired")+
  my_theme

ggsave(filename = paste0("Results/barplot_PPH4.png"),
       width = 529/72,height = 514/72,units="in",scale=1, device = "png") 

table_coloc_long <- melt(table_coloc)
table_coloc_long[, variable := gsub("posprob_coloc_", "",variable)]
table_coloc_long[, todump := ifelse(study=="meta-analysis", "blood", "liver")]
table_coloc_long[, therow := paste(hgnc, todump)]

ggplot(table_coloc_long, aes(x=clean_variable_name, y=value, fill=variable))+
  geom_bar(stat="identity", color="black")+
  facet_grid(therow ~ Category, scales = "free_x") +
  scale_fill_brewer(palette="Set1")+
  theme(axis.title.x=element_blank(),
        axis.title.y = element_blank()) +
  theme_classic() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  # scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
  # theme(strip.background = element_rect(fill = "grey90"), strip.placement = "outside", )
  theme(strip.background = element_blank(), strip.placement = "outside") +
  theme(strip.text.y = element_text(face = "bold"),
        strip.text.x = element_text(face = "bold")) +
  theme(axis.title.x = element_blank()) #+
  # my_theme

ggsave(filename = paste0("Results/barplot_coloc.png"),
       width = 529/72,height = 514/72,units="in",scale=1, device = "png") 

### forest_plot uni-cis decode and IUCPQ Biobank
unicis_forest <- res_unicis[hgnc %in% c("F2", "F11") & study %in% c("IUCPQ Biobank") & id.outcome %in% c(efficacy_outcome, safety_outcome),]
decode_forest <- res_unicis[hgnc %in% c("F2", "F11") & study %in% c("deCODE") & exposure %in% c("5316_54", "2190_55") & id.outcome %in% c(efficacy_outcome, safety_outcome),]
unicis_forest <- rbind(unicis_forest, decode_forest)

unicis_forest[, Category := ifelse(id.outcome %in%efficacy_outcome, "Efficacy outcomes", "Safety outcomes")]
unicis_forest[, exposure := gsub("SeqId_|fenland", "", exposure)]
unicis_forest[, b:=b*-1]
unicis_forest[,lci:=b-1.96*se]
unicis_forest[,uci:=b+1.96*se]
unicis_forest[, id.outcome := factor(id.outcome, levels = c(efficacy_outcome, safety_outcome))]
unicis_forest[,clean_variable_name := factor(clean_variable_name, levels = dat_balloon[order(id.outcome), unique(outcome)])]
unicis_forest <- unicis_forest[, .(Category, clean_variable_name, hgnc, study ,exposure, b, se, pval, lci, uci,
                                   posprob_coloc_PPH0, posprob_coloc_PPH1, posprob_coloc_PPH2, posprob_coloc_PPH3,posprob_coloc_PPH4), ][order(Category, clean_variable_name, hgnc, study ,exposure)]
unicis_forest[,posprob_coloc_PPH4 := round(posprob_coloc_PPH4, digits = 2)]
unicis_forest[study == "deCODE",panel:=paste0(hgnc, " blood protein levels inhibition")]
unicis_forest[study == "IUCPQ Biobank",panel:=paste0(hgnc, " liver gene levels inhibition")]
# unicis_forest[,panel:=factor(panel, levels = c("F2 liver gene levels inhibition", "F11 liver gene levels inhibition", "F2 blood protein levels inhibition", "F11 blood protein levels inhibition"))]
unicis_forest[,colour := "black"]
unicis_forest[,shape := 22]
unicis_forest[,dummy:="dummy"]
unicis_forest <- unicis_forest[order(panel, Category, clean_variable_name),]
setnames(unicis_forest, "posprob_coloc_PPH4", "PPH4")
k<-unicis_forest[study%in%"IUCPQ Biobank",]
k[,panel :=factor(panel, levels = c("F2 liver gene levels inhibition", "F11 liver gene levels inhibition"))]
k <- k[order(panel, Category, clean_variable_name),]
data <- data.table(k)
# debugonce(my_forestplotter_fancy)
p <- my_forestplotter_fancy(data = data,
                       effect.name = "Effect (95% CI)",
                       col.header = "Category",
                       col.header.heading = " ",
                       col.below.header = "clean_variable_name",
                       col.right = "PPH4",
                       col.right.heading = "PPH4",
                       exponentiate = TRUE)

ggsave(p, filename = paste0("Results/forestplot_IUCPQBiobank.png"),
       width = 1150/72,height = 350/72,units="in",scale=1, device = "png") 


k<-unicis_forest[study%in%"deCODE"]
k[,panel :=factor(panel, levels = c("F2 blood protein levels inhibition", "F11 blood protein levels inhibition"))]
k <- k[order(panel, Category, clean_variable_name),]
data <- data.table(k)
p<- my_forestplotter_fancy(data = data,
                       effect.name = "Effect (95% CI)",
                       col.header = "Category",
                       col.below.header = "clean_variable_name",
                       col.header.heading = " ",
                       col.right = "PPH4", 
                       col.right.heading = "PPH4",
                       exponentiate = TRUE)

ggsave(p, filename = paste0("Results/forestplot_unicis_decode.png"),
       width = 1150/72,height = 350/72,units="in",scale=1, device = "png") 

######ASHG presentation#######
data<- k[panel == "F2 blood protein levels inhibition", ]
data[, panel := factor(panel, levels = "F2 blood protein levels inhibition")]
p<- my_forestplotter_fancy(data = data,
                           effect.name = "Effect (95% CI)",
                           col.header = "Category",
                           col.below.header = "clean_variable_name",
                           col.header.heading = " ",
                           col.right = "PPH4", 
                           col.right.heading = "PPH4",
                           exponentiate = TRUE)
ggsave(p, filename = paste0("Results/Draft_plot/forestplot_F2_decode_ASHG.png"),
       width = 650/72,height = 250/72,units="in",scale=1, device = "png", dpi = 300) 

data<- k[panel == "F11 blood protein levels inhibition", ]
data[, panel := factor(panel, levels = "F11 blood protein levels inhibition")]
p<- my_forestplotter_fancy(data = data,
                           effect.name = "Effect (95% CI)",
                           col.header = "Category",
                           col.below.header = "clean_variable_name",
                           col.header.heading = " ",
                           col.right = "PPH4", 
                           col.right.heading = "PPH4",
                           exponentiate = TRUE)
ggsave(p, filename = paste0("Results/Draft_plot/forestplot_F11_decode_ASHG.png"),
       width = 650/72,height = 250/72,units="in",scale=1, device = "png", dpi = 300) 

####IUCPQ presentation######
data<- k[panel == "F2 blood protein levels inhibition" & clean_variable_name %in% c("Venous thromboembolism", "Cardioembolic stroke", "Bleeding", "GI-bleeding"), ]
data[, panel := factor(panel, levels = "F2 blood protein levels inhibition")]
p<- my_forestplotter_fancy(data = data,
                           effect.name = "OR (95% CI)",
                           col.header = "Category",
                           col.below.header = "clean_variable_name",
                           col.header.heading = " ",
                           # col.right = "PPH4", 
                           # col.right.heading = "PPH4",
                           exponentiate = TRUE)
ggsave(p, filename = paste0("Results/Draft_plot/forestplot_F2_decode_IUCPQ.png"),
       width = 520/72,height = 130/72, units="in",scale=1, device = "png", dpi = 300)

data<- k[panel == "F11 blood protein levels inhibition" & clean_variable_name %in% c("Venous thromboembolism", "Cardioembolic stroke", "Bleeding", "GI-bleeding"), ]
data[, panel := factor(panel, levels = "F11 blood protein levels inhibition")]
p<- my_forestplotter_fancy(data = data,
                           effect.name = "OR (95% CI)",
                           col.header = "Category",
                           col.below.header = "clean_variable_name",
                           col.header.heading = " ",
                           # col.right = "PPH4", 
                           # col.right.heading = "PPH4",
                           exponentiate = TRUE)
ggsave(p, filename = paste0("Results/Draft_plot/forestplot_F11_decode_IUCPQ.png"),
       width = 520/72,height = 130/72,units="in",scale=1, device = "png", dpi = 300) 
###### multicis decode cohorts #######
multi_forest <- res_multicis[hgnc %in% c("F2", "F11") &
                               study%in% c("deCODE") & exposure %in% c("5316_54", "2190_55") &
                               id.outcome %in% c(efficacy_outcome, safety_outcome), ]
multi_forest[, Category := ifelse(id.outcome %in% efficacy_outcome, "Efficacy outcomes", "Safety outcomes")]
multi_forest[, exposure := gsub("SeqId_|fenland", "", exposure)]
multi_forest[, b:=b*-1]
multi_forest[,lci:=b-1.96*se]
multi_forest[,uci:=b+1.96*se]
multi_forest[, id.outcome := factor(id.outcome, levels = c(efficacy_outcome, safety_outcome))]
multi_forest[,clean_variable_name := factor(clean_variable_name, levels = dat_balloon[order(id.outcome), unique(outcome)])]

multi_forest <- multi_forest[, .(Category, clean_variable_name, hgnc, study ,exposure, b, se, pval, lci, uci,
                                 nsnp.multi_cis), ][order(Category, clean_variable_name, hgnc, study ,exposure)]

multi_forest <- multi_forest[order(Category, clean_variable_name, hgnc, study ,exposure),]
multi_forest[,panel:=paste0(hgnc, " blood protein levels inhibition")]
multi_forest[,panel:=factor(panel, levels = c("F2 blood protein levels inhibition", "F11 blood protein levels inhibition"))]
multi_forest[,colour := "black"]
multi_forest[,shape := 22]
multi_forest[,dummy:="dummy"]

p <- my_forestplotter_fancy(data = multi_forest,
                       effect.name = "Effect (95% CI)",
                       col.header = "Category",
                       col.below.header = "clean_variable_name", 
                       col.header.heading = " ",
                       col.right = NULL, 
                       exponentiate = TRUE)

ggsave(p, filename = paste0("Results/forestplot_multicis_decode.png"),
       width = 1150/72,height = 350/72,units="in",scale=1, device = "png")
###### multicis multiple cohorts #######
res_multicis <- fread("Data/Modified/res_multicis_meta.txt")

multi_forest <- res_multicis[hgnc %in% c("F2", "F11") & id.outcome %in% c("dis-14-8", "dis-14-7", "dis-1-1", "dis-15-157") &
                               study%in% c("deCODE", "FENLAND", "ARIC"), ]
multi_forest[, Category := ifelse(id.outcome %in% efficacy_outcome, "Efficacy outcomes", "Safety outcomes")]
multi_forest[, exposure := gsub("SeqId_|fenland", "", exposure)]
multi_forest[, b:=b*-1]
multi_forest[,lci:=b-1.96*se]
multi_forest[,uci:=b+1.96*se]

multi_forest <- multi_forest[, .(Category, clean_variable_name, hgnc, study ,exposure, b, pval, lci, uci,
                                 nsnp.multi_cis), ][order(Category, clean_variable_name, hgnc, study ,exposure)]


dt<-data.table(Category = c(rep("Efficacy outcomes",  3), "Safety outcomes"), 
               clean_variable_name = multi_forest$clean_variable_name %>% unique,
               hgnc = "F11",
               study = "deCODE",
               exposure = "")
multi_forest <- rbind(multi_forest, dt, fill = TRUE)
multi_forest <- multi_forest[order(Category, clean_variable_name, hgnc, study ,exposure),]
multi_forest[,panel:=paste0(hgnc, " blood protein levels inhibition")]
multi_forest[,panel:=factor(panel, levels = c("F2 blood protein levels inhibition", "F11 blood protein levels inhibition"))]
multi_forest[,colour := "black"]
multi_forest[,shape := 22]
multi_forest[,dummy:="dummy"]
my_forest_fancy(
  data= multi_forest,
  col.right = c("pval"),
  col.right.heading = list(c("Effect (95% CI)", "P value")),
  exponentiate = TRUE,
  xlab="",
  col.left = c("exposure", "nsnp.multi_cis"),
  col.left.heading=c("Aptamer", " nsnp"),
  col.toformat = NULL,
  col.header = c(heading1="clean_variable_name", heading2="study", heading3="dummy"),
  yesleftcol = FALSE,
  blankrows = c(0,1,0,0))

ggsave(filename = paste0("Results/forestplot_multicis_multiplecohort.png"),
       width = 900/72,height = 350/72,units="in",scale=1, device = "png")

###### unicis multiple cohorts
uni_forest <- res_unicis[hgnc %in% c("F2", "F11") & id.outcome %in% c("dis-14-8", "dis-14-7", "dis-1-1", "dis-15-157") &
                           study%in% c("deCODE", "FENLAND", "ARIC"), ]
uni_forest[, Category := ifelse(id.outcome %in% efficacy_outcome, "Efficacy outcomes", "Safety outcomes")]
uni_forest[, exposure := gsub("SeqId_|fenland", "", exposure)]
uni_forest[, b:=b*-1]
uni_forest[,lci:=b-1.96*se]
uni_forest[,uci:=b+1.96*se]
setnames(uni_forest, "posprob_coloc_PPH4", "PPH4")
dt<-data.table(Category = "Efficacy outcomes", clean_variable_name = "Cardioembolic stroke",
               study = rep("ARIC", 2) , hgnc = "F2", exposure = c("4157_2", "5316_54"))
uni_forest <- rbind(uni_forest, dt, fill = TRUE)
uni_forest <- uni_forest[, .(Category, clean_variable_name, hgnc, study ,exposure, b, pval, lci, uci,
                             PPH4), ][order(Category, clean_variable_name, hgnc, study ,exposure)]

k1<-uni_forest[,.(hgnc, exposure)] %>% distinct
k1[, aptamerN := paste0("aptamer", 1:.N), by = "hgnc"]
uni_forest <- merge(uni_forest, k1, by = c("hgnc", "exposure"))
uni_forest <- dcast(uni_forest, Category + clean_variable_name + study + hgnc  ~ aptamerN , value.var = c("b","pval", "lci", "uci", "PPH4", "exposure"))
uni_forest <- melt(uni_forest, id.vars = c("Category", "clean_variable_name", "study", "hgnc"))
uni_forest <- separate(uni_forest, col = "variable", into = c("colnom", "todump"), sep = "_")
uni_forest <- dcast(uni_forest, Category + clean_variable_name + study  + hgnc + todump  ~  colnom, value.var = "value")

num_col<-c("b","pval", "lci", "uci")
uni_forest[, (num_col) := lapply(.SD, as.numeric),.SDcols = num_col ]
uni_forest[,PPH4:=round(as.numeric(PPH4), digits = 2)]
uni_forest <- uni_forest[order(Category, clean_variable_name, hgnc, study ,exposure),]
uni_forest[,panel:=paste0(hgnc, " blood protein levels inhibition")]
uni_forest[,panel:=factor(panel, levels = c("F2 blood protein levels inhibition", "F11 blood protein levels inhibition"))]
uni_forest[,colour := "black"]
uni_forest[,shape := 22]
uni_forest[,dummy:="dummy"]
my_forest_fancy(
  data= uni_forest,
  col.right = c("pval", "PPH4"),
  col.right.heading = list(c("Effect (95% CI)", "P value", "PPH4")),
  exponentiate = TRUE,
  xlab="",
  col.left = c("exposure"),
  col.left.heading=c("Aptamer"),
  col.toformat = NULL,
  col.header = c(heading1="clean_variable_name", heading2="study", heading3="dummy"),
  yesleftcol = FALSE,
  blankrows = c(0,1,0,0))

ggsave(filename = paste0("Results/forestplot_unicis_multiplecohort.png"),
       width = 900/72,height = 350/72,units="in",scale=1, device = "png")
######Efficacy Vs Safety for different mechanism of inhibition#####
k<-res_unicis[hgnc %in% c("F2", "F11") & id.outcome %in% c("dis-1-1","dis-15-157"),]
k[, exposure := gsub("SeqId_|fenland", "", exposure)]
k[, b:=b*-1]
k[,label:=paste0(study, "_", exposure)]
k<-k[!(study%in%c("ARIC", "FENLAND", "deCODE")),]
k[, mechanism := ifelse(study %in% c("ARIC", "FENLAND", "deCODE", "meta-analysis"), 
                        "Blood protein levels", 
                        "Hepatic gene expression")]

idtosel <- c("dis-1-1"="x", "dis-15-157"="y")
k <- k[id.outcome%in%names(idtosel), ]
k[,outcome2 := idtosel[id.outcome]]
wide <- dcast(k,mechanism + label + hgnc ~ outcome2, value.var = c("b", "se"))
wide[,mechanism := as.factor(mechanism)]
ggplot(wide, aes(x=b_x, y=b_y, color = mechanism)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = b_y - se_y, ymax = b_y + se_y),
                         colour = "grey", width = 0) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = b_x -se_x, xmax = b_x + se_x),
                          colour = "grey", height = 0)+
  facet_grid(rows = vars(hgnc)) +
  ggplot2::geom_point() +
  ggrepel::geom_text_repel(label=wide$label, check_overlap = T) +
  # scale_color_brewer("Paired") +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  ggplot2::labs(x = paste("VTE risk (log(OR)) via 1 SD decrease in coagulation factor"),
                y = paste("Bleeding risk (log(OR) via 1 SD decrease in coagulation factor")) +
  my_theme +
  geom_hline(aes(yintercept = 0),  size = 0.5) +
  geom_vline(aes(xintercept = 0),  size = 0.5) 


ggsave(filename = paste0("Results/efficacy_safety_drugmechanism.png"),
       width = 500/72,height = 500/72,units="in",scale=1, device = "png")

#######Manhattan plot
source("/mnt/sda/gagelo01/Projects/Brain_pQTL/Analysis/CMplot_eloi.R")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
# tomanhatthan<- res_map[ , .SD[which.max(posprob_coloc_PPH4)] , by = "hgnc_symbol", .SDcols = c("posprob_colocH4.SNP")]
data <- VariantAnnotation::readVcf("/mnt/sda/gagelo01/Vcffile/Server_vcf/dis-1-1/dis-1-1.vcf.gz") %>%
  gwasglue::gwasvcf_to_TwoSampleMR(., type = "outcome") %>% as.data.table(.)

data <- data[!(is.na(chr.outcome) | is.na(pos.outcome) | is.na(pval.outcome)),]
data <- data[,.(SNP, chr.outcome, pos.outcome, pval.outcome)]

setwd(paste0(wd, "/Results"))

CMplot_eloi(as.data.frame(data), type="p", plot.type="m", LOG10=TRUE, threshold=5e-08, amplify=FALSE, memo="", dpi=300, verbose=TRUE, width=14,
            height=6, col=c("grey70", "grey90"), threshold.lwd=2, cex=0.4, 
            file.output = TRUE, arrow_col_eloi = "black", highlight.text.cex=1, padding = 3, file=c("jpg"))

setwd(paste0(wd))

####QQplot
index<- sample(1:data[,.N], size = 2000, replace = FALSE)
plotdata <- data.table(
  observed = -log10(sort(data$pval.outcome)),
  expected = -log10(ppoints(data[,.N])))

qqplot <- ggplot(plotdata, aes(x = expected, y = observed)) +
  geom_point() +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               size = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -log"[10],"(", italic(P),")")),
       y = expression(paste("Observed -log"[10],"(", italic(P),")"))) +
  theme_minimal() 

# qqplot <- ggplot(data, ae(sample=pval.outcome)) +
#   stat_qq() + 
#   stat_qq_line() +
#   xlab(expression(Expected ~ ~-log[10](italic(p)))) +
#     ylab(expression(Observed ~ ~-log[10](italic(p)))) +
#   theme_minimal()

ggsave(filename = paste0("Results/qqplot_vte.png"), plot = qqplot,
       width = 350/72,height = 350/72,units="in",scale=1, device = "png")


message("This script finished without errors")


