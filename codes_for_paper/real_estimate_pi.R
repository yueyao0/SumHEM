###############################################################################
# Project: SumHEM
# Task: UKB validation | estimate pi with LDpred2-auto in whole-genome
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
phenoID <- gsub(x = args[grep(x = args, pattern = "phenoID=")], pattern = "phenoID=", replacement = "")
file_est <- gsub(x = args[grep(x = args, pattern = "file_est=")], pattern = "file_est=", replacement = "")
h2 <- gsub(x = args[grep(x = args, pattern = "h2=")], pattern = "h2=", replacement = "")
NCORES <- gsub(x = args[grep(x = args, pattern = "NCORES=")], pattern = "NCORES=", replacement = "")
path_out <- gsub(x = args[grep(x = args, pattern = "path_out=")], pattern = "path_out=", replacement = "")

h2 <- as.numeric(h2)
NCORES <- as.numeric(NCORES)

require(dplyr)
tryCatch(expr = {
  st = Sys.time()
  # GWAS
  cat(paste0("(",phenoID,") Loading GWAS data ...\n"))
  df_gwas <- readRDS(file_est) %>% select(rsid, chr, pos, beta, se, n_complete_samples)
  # LD info
  cat(paste0("(",phenoID,") Loading LD info ...\n"))
  df_map <- readRDS("/opt/working/projects/prj_035_LDpredRR/1000G/LD_EUR/map.rds") %>% select(rsid, chr, pos)
  # merge and impute missing SNP with beta=0,se=1
  df_gwas <- df_gwas %>% right_join(., df_map, by = c("rsid","chr","pos")) %>%
    rename(SNP = rsid, n_eff = n_complete_samples, beta_se = se) %>%
    mutate(na_SNP = is.na(beta), 
           beta = ifelse(is.na(beta),0,beta), 
           beta_se = ifelse(is.na(beta_se),1,beta_se), 
           n_eff = ifelse(is.na(n_eff),median(n_eff,na.rm = TRUE),n_eff), 
           scale = sqrt(beta^2+n_eff*beta_se^2), 
           sbeta = beta/scale)
  # LD matrix
  R_SFBM <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/validation/R_SFBM/R_gw.rds"))
  R_SFBM$backingfile <- paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/validation/R_SFBM/R_gw.sbk")
  # estimates
  cat(paste0("(",phenoID,") Estimating with LDpred2-auto in whole genome ...\n"))
  model <- bigsnpr::snp_ldpred2_auto(corr = R_SFBM, df_beta = df_gwas, h2_init = h2, 
                                     vec_p_init = bigsnpr::seq_log(1e-4,0.2,length.out=30), ncores = NCORES)
  # save
  saveRDS(model, file = paste0(path_out,"/model_gw.",phenoID,".rds"))
  ed = Sys.time()
  cat(paste0("(",phenoID,") Estimation is Done !!! Time ",difftime(ed,st,units="mins")," min \n"))
}, error = function(e){
  cat(e)
  ed = Sys.time()
  stop(paste0("(",phenoID,") Error !!! Time ",difftime(ed,st,units="mins")," min \n"))
})



















