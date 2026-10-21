###############################################################################
# Project: SumHEM
# Task: validation with UKB summary statistics
###############################################################################

##### ---------------------------- run in cmd --------------------------- #####
df_pheno <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/validation/res_ESTfemale_VALmale/df_pheno.rds"))

##### estimation (effects)
ws = 500
NCORES = 10
filename = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/real_estimate.sh")
if(file.exists(filename)) system(paste0("rm ",filename))
for (YY in 1:nrow(df_pheno)) {
  cm <- paste0("Rscript /opt/working/projects/prj_035_LDpredRR/UKBB/real_estimate.R", 
               " phenoID=",df_pheno$phenotype[YY], 
               " file_est=",paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/formatted_gwas/df_formatted_gwas.female.",df_pheno$phenotype[YY],".rds"), 
               " h2=",df_pheno$h2_ldsc[YY],
               " ws=",ws, 
               " NCORES=",NCORES, 
               " path_out=","/opt/working/projects/prj_035_LDpredRR/UKBB/validation/res_ESTfemale_VALmale/res_est",
               "\n")
  cat(cm, file = filename, append = TRUE)
}

##### evaluation
filename = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/real_evaluate.sh")
if(file.exists(filename)) system(paste0("rm ",filename))
for (YY in 1:nrow(df_pheno)) {
  cm <- paste0("Rscript /opt/working/projects/prj_035_LDpredRR/UKBB/real_evaluate.R", 
               " phenoID=",df_pheno$phenotype[YY], 
               " file_est=",paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/validation/res_ESTfemale_VALmale/res_est/df_est.",df_pheno$phenotype[YY],".rds"), 
               " file_val=",paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/formatted_gwas/df_formatted_gwas.male.",df_pheno$phenotype[YY],".rds"), 
               " path_out=","/opt/working/projects/prj_035_LDpredRR/UKBB/validation/res_ESTfemale_VALmale/res_val",
               "\n")
  cat(cm, file = filename, append = TRUE)
}

##### estimation (pi)
filename = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/real_estimate_pi.sh")
if(file.exists(filename)) system(paste0("rm ",filename))
for (YY in 1:nrow(df_pheno)) {
  cm <- paste0("Rscript /opt/working/projects/prj_035_LDpredRR/UKBB/real_estimate_pi.R", 
               " phenoID=",df_pheno$phenotype[YY], 
               " file_est=",paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/formatted_gwas/df_formatted_gwas.female.",df_pheno$phenotype[YY],".rds"), 
               " h2=",df_pheno$h2_ldsc[YY],
               " NCORES=",1, 
               " path_out=","/opt/working/projects/prj_035_LDpredRR/UKBB/validation/res_ESTfemale_VALmale/res_est_gw",
               "\n")
  cat(cm, file = filename, append = TRUE)
}

##### summarise
require(dplyr)
require(ggplot2)
# data
{
  df_pheno <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/validation/res_ESTfemale_VALmale/df_pheno.rds"))
  df_pi <- NULL
  for (YY in 1:nrow(df_pheno)) {
    model <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/validation/res_ESTfemale_VALmale/res_est_gw/model_gw.",df_pheno$phenotype[YY],".rds"))
    df_pi <-  rbind(df_pi, data.frame(phenotype = df_pheno$phenotype[YY], pi_est = model[[1]]$p_est))
  }
  df_eval <- NULL
  for (YY in 1:nrow(df_pheno)) {
    df_eval <- rbind(df_eval, readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/validation/res_ESTfemale_VALmale/res_val/df_eval.",df_pheno$phenotype[YY],".rds")))
  }
  df_eval <- df_eval %>% left_join(., df_pi, by=c("phenotype")) %>% left_join(., df_pheno, by=c("phenotype"))
  
}
saveRDS(df_eval, file = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/validation/res_ESTfemale_VALmale/df_eval.rds"))














