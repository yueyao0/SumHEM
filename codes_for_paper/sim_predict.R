###############################################################################
# Project: SumHEM
# Task: simulation | predicted values
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
h2 <- gsub(x = args[grep(x = args, pattern = "h2=")], pattern = "h2=", replacement = "")
pi <- gsub(x = args[grep(x = args, pattern = "pi=")], pattern = "pi=", replacement = "")
h2input <- gsub(x = args[grep(x = args, pattern = "h2input=")], pattern = "h2input=", replacement = "")
ws <- gsub(x = args[grep(x = args, pattern = "ws=")], pattern = "ws=", replacement = "")
simID <- gsub(x = args[grep(x = args, pattern = "simID=")], pattern = "simID=", replacement = "")
workpath <- gsub(x = args[grep(x = args, pattern = "workpath=")], pattern = "workpath=", replacement = "")

##### dependencies
require(doParallel)
require(dplyr)
require(Matrix)
ind_train = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_train.rds"))
ind_test = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_test.rds"))
nSNP = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_map/nSNP.rds"))
N_train = length(ind_train)
N_test = length(ind_test)
M = sum(nSNP)

##### parameters
h2 <- as.numeric(h2)
pi <- as.numeric(pi)
if (length(h2input)!=0) {
  h2input <- as.numeric(h2input)
} else {
  h2input <- h2
}
ws <- as.numeric(ws)
simID <- as.numeric(simID)

##### create direction
{
  if (h2input==h2 & ws==500) {
    if (!fs::dir_exists(paste0(workpath,"/res_pred"))) fs::dir_create(paste0(workpath,"/res_pred"))
  } else if (h2input!=h2 & ws==500) {
    if (!fs::dir_exists(paste0(workpath,"/res_pred_h2input",h2input))) fs::dir_create(paste0(workpath,"/res_pred_h2input",h2input))
  } else if (h2input==h2 & ws!=500) {
    if (!fs::dir_exists(paste0(workpath,"/res_pred_ws",ws))) fs::dir_create(paste0(workpath,"/res_pred_ws",ws))
  }
}

##### load the data
{
  if (h2input==h2 & ws==500) {
    df_est <- readRDS(paste0(workpath,"/res_est/df_est_simID",simID,".rds"))
  } else if (h2input!=h2 & ws==500) {
    df_est <- readRDS(paste0(workpath,"/res_est_h2input",h2input,"/df_est_simID",simID,".rds"))
  } else if (h2input==h2 & ws!=500) {
    df_est <- readRDS(paste0(workpath,"/res_est_ws",ws,"/df_est_simID",simID,".rds"))
  }
}

##### prediction
st_ID = Sys.time()
{
  cat(paste0("========== prediction h2=",h2," pi=",round(pi,4)," h2input=",h2input," ws=",ws," simID=",simID," ==========","\n"))
  ls_y_pred <- list(simID = simID, 
                    y_true = numeric(N_train+N_test), y_gwas = numeric(N_train+N_test), 
                    y_inf_LDpred = numeric(N_train+N_test), y_auto_LDpred = numeric(N_train+N_test), 
                    y_1_RRpred = numeric(N_train+N_test), y_2_RRpred = numeric(N_train+N_test))
  for (CHR in 1:22) {
    cat(paste0("prediction CHR ", CHR, "/22","\t"))
    st_chr = Sys.time()
    obj.bigSNP <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_imputed/chr_",CHR,".rds"))
    obj.bigSNP$genotypes$backingfile <- paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_imputed/chr_",CHR,".bk")
    obj.bigSNP$genotypes$code256 <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_imputed_codes/chr_",CHR,".rds"))$code256
    sG_train <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_scaled/sG_train_chr", CHR, ".rds"))
    sG_test <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_scaled/sG_test_chr", CHR, ".rds"))
    df_est_chr <- df_est %>% filter(chr == CHR)
    # true effects
    ls_y_pred$y_true[ind_train] <- ls_y_pred$y_true[ind_train] + 
      bigstatsr::big_prodVec(X = obj.bigSNP$genotypes, ind.row = ind_train, y.col = df_est_chr$b_true, 
                             center = sG_train$center, scale = sG_train$scale, ncores = 1)
    ls_y_pred$y_true[ind_test] <- ls_y_pred$y_true[ind_test] + 
      bigstatsr::big_prodVec(X = obj.bigSNP$genotypes, ind.row = ind_test, y.col = df_est_chr$b_true, 
                             center = sG_test$center, scale = sG_test$scale, ncores = 1)
    # gwas
    ls_y_pred$y_gwas[ind_train] <- ls_y_pred$y_gwas[ind_train] + 
      bigstatsr::big_prodVec(X = obj.bigSNP$genotypes, ind.row = ind_train, y.col = df_est_chr$sbeta, 
                             center = sG_train$center, scale = sG_train$scale, ncores = 1)
    ls_y_pred$y_gwas[ind_test] <- ls_y_pred$y_gwas[ind_test] + 
      bigstatsr::big_prodVec(X = obj.bigSNP$genotypes, ind.row = ind_test, y.col = df_est_chr$sbeta, 
                             center = sG_test$center, scale = sG_test$scale, ncores = 1)
    # LDpred-inf
    ls_y_pred$y_inf_LDpred[ind_train] <- ls_y_pred$y_inf_LDpred[ind_train] + 
      bigstatsr::big_prodVec(X = obj.bigSNP$genotypes, ind.row = ind_train, y.col = df_est_chr$b_inf_LDpred, 
                             center = sG_train$center, scale = sG_train$scale, ncores = 1)
    ls_y_pred$y_inf_LDpred[ind_test] <- ls_y_pred$y_inf_LDpred[ind_test] + 
      bigstatsr::big_prodVec(X = obj.bigSNP$genotypes, ind.row = ind_test, y.col = df_est_chr$b_inf_LDpred, 
                             center = sG_test$center, scale = sG_test$scale, ncores = 1)
    # LDpred-auto
    ls_y_pred$y_auto_LDpred[ind_train] <- ls_y_pred$y_auto_LDpred[ind_train] + 
      bigstatsr::big_prodVec(X = obj.bigSNP$genotypes, ind.row = ind_train, y.col = df_est_chr$b_auto_LDpred, 
                             center = sG_train$center, scale = sG_train$scale, ncores = 1)
    ls_y_pred$y_auto_LDpred[ind_test] <- ls_y_pred$y_auto_LDpred[ind_test] + 
      bigstatsr::big_prodVec(X = obj.bigSNP$genotypes, ind.row = ind_test, y.col = df_est_chr$b_auto_LDpred, 
                             center = sG_test$center, scale = sG_test$scale, ncores = 1)
    # SumRR
    ls_y_pred$y_SumRR[ind_train] <- ls_y_pred$y_SumRR[ind_train] + 
      bigstatsr::big_prodVec(X = obj.bigSNP$genotypes, ind.row = ind_train, y.col = df_est_chr$b_SumRR, 
                             center = sG_train$center, scale = sG_train$scale, ncores = 1)
    ls_y_pred$y_SumRR[ind_test] <- ls_y_pred$y_SumRR[ind_test] + 
      bigstatsr::big_prodVec(X = obj.bigSNP$genotypes, ind.row = ind_test, y.col = df_est_chr$b_SumRR, 
                             center = sG_test$center, scale = sG_test$scale, ncores = 1)
    # SumHEM
    ls_y_pred$y_SumHEM[ind_train] <- ls_y_pred$y_SumHEM[ind_train] + 
      bigstatsr::big_prodVec(X = obj.bigSNP$genotypes, ind.row = ind_train, y.col = df_est_chr$b_SumHEM, 
                             center = sG_train$center, scale = sG_train$scale, ncores = 1)
    ls_y_pred$y_SumHEM[ind_test] <- ls_y_pred$y_SumHEM[ind_test] + 
      bigstatsr::big_prodVec(X = obj.bigSNP$genotypes, ind.row = ind_test, y.col = df_est_chr$b_SumHEM, 
                             center = sG_test$center, scale = sG_test$scale, ncores = 1)
    rm(obj.bigSNP, sG_train, sG_test, df_est_chr)
    ed_chr = Sys.time()
    cat(paste0("Done !!! Time ",difftime(ed_chr,st_chr,units="mins")," min \n"))
  }
  # save
  if (h2input==h2 & ws==500) {
    saveRDS(ls_y_pred, file = paste0(workpath,"/res_pred/ls_y_pred_simID",simID,".rds"))
  } else if (h2input!=h2 & ws==500) {
    saveRDS(ls_y_pred, file = paste0(workpath,"/res_pred_h2input",h2input,"/ls_y_pred_simID",simID,".rds"))
  } else if (h2input==h2 & ws!=500) {
    saveRDS(ls_y_pred, file = paste0(workpath,"/res_pred_ws",ws,"/ls_y_pred_simID",simID,".rds"))
  }
  rm(df_est,ls_y_pred)
}
ed_ID = Sys.time()
cat(paste0("Prediction is done !!! Time ",difftime(ed_ID,st_ID,units="mins")," min \n"))

