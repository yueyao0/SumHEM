###############################################################################
# Project: SumHEM
# Task: UKB validation | evaluate r2 and h2 in summary-level
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
phenoID <- gsub(x = args[grep(x = args, pattern = "phenoID=")], pattern = "phenoID=", replacement = "")
file_est <- gsub(x = args[grep(x = args, pattern = "file_est=")], pattern = "file_est=", replacement = "")
file_val <- gsub(x = args[grep(x = args, pattern = "file_val=")], pattern = "file_val=", replacement = "")
path_out <- gsub(x = args[grep(x = args, pattern = "path_out=")], pattern = "path_out=", replacement = "")

require(dplyr)
tryCatch(expr = {
  st = Sys.time()
  # data
  cat(paste0("(",phenoID,") Loading estimation results ...\n"))
  df_est <- readRDS(file_est)
  cat(paste0("(",phenoID,") Loading validation GWAS ...\n"))
  df_val <- readRDS(file_val) %>%
    select(rsid, chr, pos, beta, se, n_complete_samples) %>%
    rename(SNP = rsid, n_eff = n_complete_samples, beta_se = se) %>%
    mutate(na_SNP = is.na(beta), 
           beta = ifelse(is.na(beta),0,beta), 
           beta_se = ifelse(is.na(beta_se),1,beta_se), 
           n_eff = ifelse(is.na(n_eff),median(n_eff,na.rm=TRUE),n_eff), 
           scale = sqrt(beta^2+n_eff*beta_se^2), 
           sbeta = beta/scale)
  # r2
  cat(paste0("(",phenoID,") Evaluating ...\n"))
  df_eval <- data.frame(phenotype = phenoID, 
                        h2_inf_LDpred = as.numeric(df_est$b_inf_LDpred %*% df_est$Rb_inf_LDpred), 
                        h2_auto_LDpred = as.numeric(df_est$b_auto_LDpred %*% df_est$Rb_auto_LDpred), 
                        h2_SumRR = as.numeric(df_est$b_SumRR %*% df_est$Rb_SumRR), 
                        h2_SumHEM = as.numeric(df_est$b_SumHEM %*% df_est$Rb_SumHEM), 
                        r2_est_inf_LDpred = as.numeric(df_est$b_inf_LDpred %*% df_est$sbeta), 
                        r2_est_auto_LDpred = as.numeric(df_est$b_auto_LDpred %*% df_est$sbeta), 
                        r2_est_SumRR = as.numeric(df_est$b_SumRR %*% df_est$sbeta), 
                        r2_est_SumHEM = as.numeric(df_est$b_SumHEM %*% df_est$sbeta), 
                        r2_val_inf_LDpred = as.numeric(df_est$b_inf_LDpred %*% df_val$sbeta)^2/as.numeric(df_est$b_inf_LDpred %*% df_est$sbeta), 
                        r2_val_auto_LDpred = as.numeric(df_est$b_auto_LDpred %*% df_val$sbeta)^2/as.numeric(df_est$b_auto_LDpred %*% df_est$sbeta), 
                        r2_val_SumRR = as.numeric(df_est$b_SumRR %*% df_val$sbeta)^2/as.numeric(df_est$b_SumRR %*% df_est$sbeta), 
                        r2_val_SumHEM = as.numeric(df_est$b_SumHEM %*% df_val$sbeta)^2/as.numeric(df_est$b_SumHEM %*% df_est$sbeta))
  # save
  saveRDS(df_eval, file = paste0(path_out,"/df_eval.",phenoID,".rds"))
  rm(df_est,df_val)
  ed = Sys.time()
  cat(paste0("(",phenoID,") Evaluation is Done !!! Time ",difftime(ed,st,units="mins")," min \n"))
}, error = function(e){
  cat(e)
  ed = Sys.time()
  stop(paste0("(",phenoID,") Error !!! Time ",difftime(ed,st,units="mins")," min \n"))
})














