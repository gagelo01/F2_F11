wald_quick <- function(inst_all_sign_clump_arg, all_outcome_arg, id_exposure_name, id_outcome_name = ".*") {

  message(paste0("###### Analysing ", id_exposure_name, " #########"))

  harm_univariate <- TwoSampleMR::harmonise_data(exposure_dat = inst_all_sign_clump_arg[grep(paste0("^" ,id_exposure_name, "$"), id.exposure),],
                                                 outcome_dat = all_outcome_arg[ grep(paste0("^" ,id_outcome_name, "$"), id.outcome),],
                                                 action = 1)

  if (dim(harm_univariate)[1] == 0) {
    return(cbind(inst_all_sign_clump_arg[grep(paste0("^" ,id_exposure_name, "$"), id.exposure), .(id.exposure,exposure)],
                 dplyr::distinct(all_outcome_arg[grep(paste0("^" ,id_outcome_name, "$"), id.outcome), .(id.outcome,outcome)])))
  }

  harm_univariate <- TwoSampleMR::steiger_filtering(harm_univariate)
  data.table::setDT(harm_univariate)
  harm_univariate <- harm_univariate[mr_keep==TRUE,]
  res <- TwoSampleMR::mr(harm_univariate, method_list = c("mr_wald_ratio"))
  res$nsnp <- harm_univariate$SNP
  res$type_of_test <- "Primary analysis"
  data.table::setDT(res)
  res$lci<-(res$b)-(res$se*1.96)
  res$uci <- (res$b)+(res$se*1.96)
  harm_univariate <- TwoSampleMR::add_rsq(harm_univariate) %>% data.table::as.data.table(.)
  res$fstat<-GagnonMR::fstat_fromdat(split(harm_univariate, 1:nrow(harm_univariate)))
  res$rsq <- harm_univariate$rsq.exposure
  res <- merge(res, harm_univariate[,.(id.exposure,id.outcome,steiger_dir,steiger_pval)], by = c("id.exposure", "id.outcome"))
  data.table::setDT(res)
  return(res)
}
