###############################################################################
# Project: SumHEM
# Task: simulation with UKB genotype
###############################################################################

##### split dataset
N = 336000 # 200000 for train, 100000 for test
ind_train <- sort(sample(1:N, size = round(N/2)))
ind_test <- sort(setdiff(1:N, ind_train))
saveRDS(ind_train, file = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_train.rds"))
saveRDS(ind_test, file = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_test.rds"))

##### genotype: missing_rate<0.02, MAF>0.01, Hardy-Weinberg equilibrium p-value<1e-5, impute NA via bigsnpr::snp_fastImputeSimple (random)
# ./plink --bfile ./ukb_unrelated_british_n336000/ukb_chr1_n336000 --geno 0.02 --maf 0.01 --hwe 1e-5 --make-bed -out ./0.01/qc0.01/0.01chr_1

##### R matrix (train)
{
  # plink remove SNP with missing rate
  for (CHR in 1:22) {
    cmd <- paste0("/opt/working/apps/plink",
                  " --bfile ","/opt/storage/ukb_unrelated_british_n336000/ukb_chr",CHR,"_n336000",
                  " --geno 0.02 --maf 0.01 --hwe 1e-5",
                  " --make-bed",
                  " --out ","/prj_035_LDpredRR/UKBB/simulation/G_filter/G_filter_chr",CHR,"_n336000",
                  "\n")
    if (CHR==1) cat(cmd, file = "/Users/yyds/Documents/Projects/PRS/RRpred/UKBB/plink_G.txt", append = FALSE)
    else cat(cmd, file = "/Users/yyds/Documents/Projects/PRS/RRpred/UKBB/plink_G.txt", append = TRUE)
  }
}

require(dplyr)
ind_train = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_train.rds"))
ind_test = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_test.rds"))
st_res = Sys.time()
nSNP <- numeric(22)
ls_map <- vector(mode = "list", length = 22)
for (CHR in 1:22) {
  print(paste0("===== preparation scale and R matrix CHR ", CHR, "/22 ====="))
  st_chr = Sys.time()
  print("loading")
  obj.bigSNP <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_imputed/chr_",CHR,".rds"))
  obj.bigSNP$genotypes$backingfile <- paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_imputed/chr_",CHR,".bk")
  obj.bigSNP$genotypes$code256 <- readRDS(paste0("/opt/working/projects/prj_038_HEM/workspace/0.01/imputed_backing/chr_",CHR,".rds"))$code256
  print("beta information")
  ls_map[[CHR]] <- obj.bigSNP$map
  saveRDS(obj.bigSNP$map, file = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_map/map_chr", CHR, ".rds"))
  print("scaling")
  nSNP[CHR] <- obj.bigSNP$genotypes$ncol
  sG_train <- bigstatsr::big_scale(center = TRUE, scale = TRUE)(X = obj.bigSNP$genotypes, ind.row = ind_train)
  sG_train$scale <- ifelse(sG_train$scale==0, 1e-10, sG_train$scale)
  sG_test  <- bigstatsr::big_scale(center = TRUE, scale = TRUE)(X = obj.bigSNP$genotypes, ind.row = ind_test)
  sG_test$scale <- ifelse(sG_test$scale==0, 1e-10, sG_test$scale)
  saveRDS(nSNP, file = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_map/nSNP.rds"))
  saveRDS(sG_train, file = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_scaled/sG_train_chr", CHR, ".rds"))
  saveRDS(sG_test, file = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_scaled/sG_test_chr", CHR, ".rds"))
  print("computing R matrix")
  R <- bigsnpr::snp_cor(Gna = obj.bigSNP$genotypes, ind.row = ind_train, size = 500, ncores = 10)
  R@x[is.na(R@x)] <- 1e-10
  saveRDS(R, file = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/R_dsCM/R_chr", CHR, ".rds"))
  rm(obj.bigSNP, sG_train, sG_test, R)
  ed_chr = Sys.time()
  print(ed_chr-st_chr)
}
df_map <- data.table::rbindlist(ls_map) %>%
  rename(chr = chromosome, SNP = marker.ID, pos = physical.pos, A1 = allele1, A2 = allele2) %>%
  select(SNP, chr, pos, A1, A2) %>%
  arrange(chr, pos)
saveRDS(as.data.frame(df_map), file = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_map/df_map.rds"))
ed_res = Sys.time()
print(ed_res-st_res)

##### R matrix (test)
ind_test = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_test.rds"))
st_res = Sys.time()
for (CHR in 1:22) {
  print(paste0("===== preparation scale and R matrix CHR ", CHR, "/22 ====="))
  st_chr = Sys.time()
  print("loading")
  obj.bigSNP <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_imputed/chr_",CHR,".rds"))
  obj.bigSNP$genotypes$backingfile <- paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_imputed/chr_",CHR,".bk")
  obj.bigSNP$genotypes$code256 <- readRDS(paste0("/opt/working/projects/prj_038_HEM/workspace/0.01/imputed_backing/chr_",CHR,".rds"))$code256
  print("computing R matrix")
  R <- bigsnpr::snp_cor(Gna = obj.bigSNP$genotypes, ind.row = ind_test, size = 500, ncores = 10)
  R@x[is.na(R@x)] <- 1e-10
  saveRDS(R, file = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/R_dsCM_test/R_chr", CHR, ".rds"))
  rm(obj.bigSNP, R)
  ed_chr = Sys.time()
  print(ed_chr-st_chr)
}
ed_res = Sys.time()
print(ed_res-st_res)

##### LD-score
require(dplyr)
st_res = Sys.time()
df_map <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_map/df_map.rds"))
LDscore <- NULL
for (CHR in 1:22) {
  print(paste0("LDscore CHR ", CHR, "/22"))
  df_map_chr <- df_map %>% filter(chr == CHR)
  R_chr <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/R_dsCM/R_chr", CHR, ".rds"))
  R_chr@x[is.na(R_chr@x)] <- 1e-10
  LDscore <- c(LDscore, Matrix::colSums(R_chr^2))
  rm(df_map_chr, R_chr)
}
df_LD <- df_map %>% mutate(LDscore = LDscore)
saveRDS(df_LD, file = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/R_dsCM/df_LD.rds"))
ed_res = Sys.time()
print(ed_res-st_res)

###############################################################################
########## ------------------------ CMD ---------------------------- ##########
###############################################################################

##### simulate & GWAS (pi)
ws = 500
Nsim = 20
for (h2 in c(0.1,0.5,0.7)) {
  for (pi in c(0.1,0.3,0.5,0.7)) {
    # h2 = 0.3 # 0.1,0.5,0.7
    # pi = 0.1 # 0.1,0.3,0.5,0.7
    workpath = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/continuous/","h2",h2,"_pi",pi)
    filename = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulate","_h2",h2,"_pi",pi,".sh")
    if(file.exists(filename)) system(paste0("rm ",filename))
    for (simID in 1:Nsim) {
      cm <- paste0("Rscript /opt/working/projects/prj_035_LDpredRR/UKBB/sim_simulate.R", 
                   " h2=",h2, " pi=",pi, " ws=",ws, " simID=",simID, 
                   " workpath=",workpath, 
                   "\n")
      cat(cm, file = filename, append = TRUE)
    }
  }
}

##### estimate (pi | h2input | ws)
ws = 500
NCORES = 10
Nsim = 20
for (h2 in c(0.1,0.5,0.7)) {
  for (pi in c(0.1,0.3,0.5,0.7)) {
    # h2 = 0.5 # 0.1,0.5,0.7
    # pi = 0.1 # 0.1,0.3,0.5,0.7
    h2input = h2
    workpath = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/continuous/","h2",h2,"_pi",pi)
    filename = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/estimate","_h2",h2,"_pi",pi,".sh")
    if(file.exists(filename)) system(paste0("rm ",filename))
    for (simID in 1:Nsim) {
      cm <- paste0("Rscript /opt/working/projects/prj_035_LDpredRR/UKBB/sim_estimate.R", 
                   " h2=",h2, " pi=",pi, " h2input=",h2input, " simID=",simID, 
                   " ws=",ws, " NCORES=",NCORES, 
                   " workpath=",workpath, 
                   "\n")
      cat(cm, file = filename, append = TRUE)
    }
  }
}

##### predict (pi | h2input | ws)
ws = 500
NCORES = 10
Nsim = 20
for (h2 in c(0.1,0.5,0.7)) {
  for (pi in c(0.1,0.3,0.5,0.7)) {
    # h2 = 0.5 # 0.1,0.5,0.7
    # pi = 0.1 # 0.1,0.3,0.5,0.7
    h2input = h2
    workpath = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/continuous/","h2",h2,"_pi",pi)
    filename = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/predict_h2",h2,"_pi",pi,".sh")
    if(file.exists(filename)) system(paste0("rm ",filename))
    for (simID in 1:Nsim) {
      cm <- paste0("Rscript /opt/working/projects/prj_035_LDpredRR/UKBB/sim_predict.R", 
                   " h2=",h2, " pi=",pi, " h2input=",h2input, " ws=",ws, " simID=",simID, 
                   " workpath=",workpath, 
                   "\n")
      cat(cm, file = filename, append = TRUE)
    }
  }
}

##### estimate (LDpred2-auto without h2 updating)
ws = 500
NCORES = 10
Nsim = 20
h2 = 0.5
for (pi in c(0.1,0.3,0.5,0.7)) {
  # pi = 0.1
  h2input = h2
  workpath = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/continuous/","h2",h2,"_pi",pi)
  filename = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/sim_estimate_LDpred_no_h2update.sh")
  if(file.exists(filename)) system(paste0("rm ",filename))
  for (simID in 1:5) {
    cm <- paste0("Rscript /opt/working/projects/prj_035_LDpredRR/UKBB/sim_estimate_LDpred_no_h2update.R", 
                 " h2=",h2, " pi=",pi, " h2input=",h2input, " ws=",ws, " simID=",simID, 
                 " ws=",ws, " NCORES=",NCORES, 
                 " workpath=",workpath, 
                 "\n")
    cat(cm, file = filename, append = TRUE)
  }
}

###############################################################################
########## ---------------------- evaluation ----------------------- ##########
###############################################################################

##### evaluation: r2
ind_train = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_train.rds"))
ind_test = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_test.rds"))
nSNP = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_map/nSNP.rds"))
N_train = length(ind_train)
N_test = length(ind_test)
M = sum(nSNP)
require(dplyr)
h2 = 0.5
pi = 0.7
h2input = h2
ws = 500
Nsim = 20
for (h2 in c(0.1,0.5,0.7)) {
  for (pi in c(0.1,0.3,0.5,0.7)) {
    # for (h2input in c(0.3,0.4,0.5,0.6,0.7)) {
      # for (ws in c(250,500,1000)) {
      workpath = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/continuous/","h2",h2,"_pi",pi)
      Mcausal = round(M*pi)
      h2input = h2
      df_eval_r2 <- NULL
      for (simID in 1:Nsim) {
        print(paste0("===== evaluation: r2 | h2=",h2," pi=",pi," h2input=",h2input," ws=",ws," simID=",simID,"/",Nsim," ====="))
        # load the data
        {
          y_sim <- readRDS(paste0(workpath,"/res_true/y_sim_simID",simID,".rds"))
          if (h2input==h2 & ws==500) {
            ls_y_pred <- readRDS(paste0(workpath,"/res_pred/ls_y_pred_simID",simID,".rds"))
          } else if (h2input!=h2 & ws==500) {
            ls_y_pred <- readRDS(paste0(workpath,"/res_pred_h2input",h2input,"/ls_y_pred_simID",simID,".rds"))
          } else if (h2input==h2 & ws!=500) {
            ls_y_pred <- readRDS(paste0(workpath,"/res_pred_ws",ws,"/ls_y_pred_simID",simID,".rds"))
          }
          ls_y_pred$y_sim <- y_sim
        }
        # dataframe
        {
          df <- data.frame(simID = simID, Nsim = Nsim, N_train = N_train, N_test = N_test, M = M, h2 = h2, pi = pi, Mcausal = Mcausal, h2input = h2input, ws = ws, 
                           r2est_true = summary(lm(ls_y_pred$y_true[ind_train] ~ ls_y_pred$y_sim[ind_train]))$r.squared, 
                           r2est_gwas = summary(lm(ls_y_pred$y_gwas[ind_train] ~ ls_y_pred$y_sim[ind_train]))$r.squared, 
                           r2est_inf_LDpred = summary(lm(ls_y_pred$y_inf_LDpred[ind_train] ~ ls_y_pred$y_sim[ind_train]))$r.squared, 
                           r2est_auto_LDpred = summary(lm(ls_y_pred$y_auto_LDpred[ind_train] ~ ls_y_pred$y_sim[ind_train]))$r.squared, 
                           r2est_SumRR = summary(lm(ls_y_pred$y_SumRR[ind_train] ~ ls_y_pred$y_sim[ind_train]))$r.squared, 
                           r2est_SumHEM = summary(lm(ls_y_pred$y_SumHEM[ind_train] ~ ls_y_pred$y_sim[ind_train]))$r.squared, 
                           r2val_true = summary(lm(ls_y_pred$y_true[ind_test] ~ ls_y_pred$y_sim[ind_test]))$r.squared, 
                           r2val_gwas = summary(lm(ls_y_pred$y_gwas[ind_test] ~ ls_y_pred$y_sim[ind_test]))$r.squared, 
                           r2val_inf_LDpred = summary(lm(ls_y_pred$y_inf_LDpred[ind_test] ~ ls_y_pred$y_sim[ind_test]))$r.squared, 
                           r2val_auto_LDpred = summary(lm(ls_y_pred$y_auto_LDpred[ind_test] ~ ls_y_pred$y_sim[ind_test]))$r.squared, 
                           r2val_SumRR = summary(lm(ls_y_pred$y_SumRR[ind_test] ~ ls_y_pred$y_sim[ind_test]))$r.squared, 
                           r2val_SumHEM = summary(lm(ls_y_pred$y_SumHEM[ind_test] ~ ls_y_pred$y_sim[ind_test]))$r.squared)
        }
        df_eval_r2 <- rbind(df_eval_r2, df)
      }
      # save
      {
        if (h2input==h2 & ws==500) {
          saveRDS(df_eval_r2, file = paste0(workpath,"/df_eval_r2.rds"))
        } else if (h2input!=h2 & ws==500) {
          saveRDS(df_eval_r2, file = paste0(workpath,"/df_eval_r2_h2input",h2input,".rds"))
        } else if (h2input==h2 & ws!=500) {
          saveRDS(df_eval_r2, file = paste0(workpath,"/df_eval_r2_ws",ws,".rds"))
        }
      }
      # }
    # }
  }
}

##### evaluation: h2
ind_train = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_train.rds"))
ind_test = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_test.rds"))
nSNP = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_map/nSNP.rds"))
N_train = length(ind_train)
N_test = length(ind_test)
M = sum(nSNP)
require(dplyr)
h2 = 0.5
pi = 0.7
h2input = h2
ws = 500
Nsim = 20
for (h2 in c(0.1,0.5,0.7)) {
  for (pi in c(0.1,0.3,0.5,0.7)) {
    # for (h2input in c(0.3,0.4,0.5,0.6,0.7)) {
      # for (ws in c(250,500,1000)) {
      workpath = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/continuous/","h2",h2,"_pi",pi)
      Mcausal = round(M*pi)
      h2input = h2
      df_eval_h2 <- NULL
      for (simID in 1:Nsim) {
        print(paste0("===== evaluation: h2 | h2=",h2," pi=",pi," h2input=",h2input," ws=",ws," simID=",simID,"/",Nsim," ====="))
        # load the data
        {
          if (h2input==h2 & ws==500) {
            df_est <- readRDS(paste0(workpath,"/res_est/df_est_simID",simID,".rds"))
          } else if (h2input!=h2 & ws==500) {
            df_est <- readRDS(paste0(workpath,"/res_est_h2input",h2input,"/df_est_simID",simID,".rds"))
          } else if (h2input==h2 & ws!=500) {
            df_est <- readRDS(paste0(workpath,"/res_est_ws",ws,"/df_est_simID",simID,".rds"))
          }
        }
        # dataframe
        {
          df <- data.frame(simID = simID, Nsim = Nsim, N_train = N_train, N_test = N_test, M = M, h2 = h2, pi = pi, Mcausal = Mcausal, h2input = h2input, ws = ws, 
                           h2_true = as.numeric(df_est$b_true %*% df_est$Rb_true),
                           h2_gwas = as.numeric(df_est$sbeta %*% df_est$Rb_gwas),
                           h2_inf_LDpred = as.numeric(df_est$b_inf_LDpred %*% df_est$Rb_inf_LDpred), 
                           h2_auto_LDpred = as.numeric(df_est$b_auto_LDpred %*% df_est$Rb_auto_LDpred), 
                           h2_SumRR = as.numeric(df_est$b_SumRR %*% df_est$Rb_SumRR), 
                           h2_SumHEM = as.numeric(df_est$b_SumHEM %*% df_est$Rb_SumHEM))
        }
        df_eval_h2 <- rbind(df_eval_h2, df)
      }
      # save
      {
        if (h2input==h2 & ws==500) {
          saveRDS(df_eval_h2, file = paste0(workpath,"/df_eval_h2.rds"))
        } else if (h2input!=h2 & ws==500) {
          saveRDS(df_eval_h2, file = paste0(workpath,"/df_eval_h2_h2input",h2input,".rds"))
        } else if (h2input==h2 & ws!=500) {
          saveRDS(df_eval_h2, file = paste0(workpath,"/df_eval_h2_ws",ws,".rds"))
        }
      }
      # }
    # }
  }
}

##### evaluation: effect direction f1-score
ind_train = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_train.rds"))
ind_test = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_test.rds"))
nSNP = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_map/nSNP.rds"))
N_train = length(ind_train)
N_test = length(ind_test)
M = sum(nSNP)
require(dplyr)
h2 = 0.5
pi = 0.7
h2input = h2
ws = 500
Nsim = 20
con_mat <- function(b1, b2){
  mat <- matrix(NA, nrow = 3, ncol = 3) # col: b1 (ref), row: b2 (pred)
  mat[1,1] <- sum(b1<(-1e-5) & b2<(-1e-5))
  mat[2,1] <- sum(b1<(-1e-5) & b2>(-1e-5) & b2<1e-5)
  mat[3,1] <- sum(b1<(-1e-5) & b2>1e-5)
  mat[1,2] <- sum(b1>(-1e-5) & b1<1e-5 & b2<(-1e-5))
  mat[2,2] <- sum(b1>(-1e-5) & b1<1e-5 & b2>(-1e-5) & b2<1e-5)
  mat[3,2] <- sum(b1>(-1e-5) & b1<1e-5 & b2>1e-5)
  mat[1,3] <- sum(b1>1e-5 & b2<(-1e-5))
  mat[2,3] <- sum(b1>1e-5 & b2>(-1e-5) & b2<1e-5)
  mat[3,3] <- sum(b1>1e-5 & b2>1e-5)
  return(mat)
}
f1_score <- function(mat){
  # class 1
  TP_1 <- mat[1,1]
  FP_1 <- mat[1,2] + mat[1,3]
  FN_1 <- mat[2,1] + mat[3,1]
  TN_1 <- mat[2,2] + mat[2,3] + mat[3,2] + mat[3,3]
  P_1 <- TP_1/(TP_1+FP_1)
  R_1 <- TP_1/(TP_1+FN_1)
  f1_1 <- 2*P_1*R_1/(P_1+R_1)
  # class 2
  TP_2 <- mat[2,2]
  FP_2 <- mat[2,1] + mat[2,3]
  FN_2 <- mat[1,2] + mat[3,2]
  TN_2 <- mat[1,1] + mat[1,3] + mat[3,1] + mat[3,3]
  P_2 <- TP_2/(TP_2+FP_2)
  R_2 <- TP_2/(TP_2+FN_2)
  f1_2 <- 2*P_2*R_2/(P_2+R_2)
  # class 3
  TP_3 <- mat[3,3]
  FP_3 <- mat[3,1] + mat[3,2]
  FN_3 <- mat[1,3] + mat[2,3]
  TN_3 <- mat[1,1] + mat[1,2] + mat[2,1] + mat[2,2]
  P_3 <- TP_3/(TP_3+FP_3)
  R_3 <- TP_3/(TP_3+FN_3)
  f1_3 <- 2*P_3*R_3/(P_3+R_3)
  # micro
  micro_P <- (TP_1+TP_2+TP_3)/(TP_1+TP_2+TP_3+FP_1+FP_2+FP_3)
  micro_R <- (TP_1+TP_2+TP_3)/(TP_1+TP_2+TP_3+FN_1+FN_2+FN_3)
  micro_f1 <- 2*micro_P*micro_R/(micro_P+micro_R)
  # macro
  macro_f1 <- (f1_1+f1_2+f1_3)/3
  # weighted
  weighted_f1 <- (sum(mat[,1])*f1_1+sum(mat[,2])*f1_2+sum(mat[,3])*f1_3)/sum(mat)
  # summarise
  f1 <- list(micro_f1 = micro_f1, macro_f1 = macro_f1, weighted_f1 = weighted_f1)
  # f1 <- c(micro_f1, macro_f1, weighted_f1)
  # names(f1) <- c("micro_f1", "macro_f1", "weighted_f1")
  return(f1)
}
for (h2 in c(0.1,0.5,0.7)) {
  for (pi in c(0.1,0.3,0.5,0.7)) {
    # for (h2input in c(0.3,0.4,0.5,0.6,0.7)) {
      # for (ws in c(250,500,1000)) {
      workpath = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/continuous/","h2",h2,"_pi",pi)
      Mcausal = round(M*pi)
      h2input = h2
      ls_mat_b_true_true = ls_mat_b_true_gwas = ls_mat_b_true_infLDpred = ls_mat_b_true_autoLDpred = ls_mat_b_true_SumRR = ls_mat_b_true_SumHEM <- vector(mode = "list", length = Nsim)
      ls_f1_b_true_true = ls_f1_b_true_gwas = ls_f1_b_true_infLDpred = ls_f1_b_true_autoLDpred = ls_f1_b_true_SumRR = ls_f1_b_true_SumHEM <- vector(mode = "list", length = Nsim)
      df_eval_f1 <- NULL
      for (simID in 1:Nsim) {
        print(paste0("===== evaluation: effect direction f1-score | h2=",h2," pi=",pi," h2input=",h2input," ws=",ws," simID=",simID,"/",Nsim," ====="))
        # load the data
        {
          if (h2input==h2 & ws==500) {
            df_est <- readRDS(paste0(workpath,"/res_est/df_est_simID",simID,".rds"))
          } else if (h2input!=h2 & ws==500) {
            df_est <- readRDS(paste0(workpath,"/res_est_h2input",h2input,"/df_est_simID",simID,".rds"))
          } else if (h2input==h2 & ws!=500) {
            df_est <- readRDS(paste0(workpath,"/res_est_ws",ws,"/df_est_simID",simID,".rds"))
          }
        }
        # confusion matrix of b_est vs true
        {
          ls_mat_b_true_true[[simID]] <- con_mat(b1 = df_est$b_true, b2 = df_est$b_true)
          ls_mat_b_true_gwas[[simID]] <- con_mat(b1 = df_est$b_true, b2 = df_est$sbeta)
          ls_mat_b_true_infLDpred[[simID]] <- con_mat(b1 = df_est$b_true, b2 = df_est$b_inf_LDpred)
          ls_mat_b_true_autoLDpred[[simID]] <- con_mat(b1 = df_est$b_true, b2 = df_est$b_auto_LDpred)
          ls_mat_b_true_SumRR[[simID]] <- con_mat(b1 = df_est$b_true, b2 = df_est$b_SumRR)
          ls_mat_b_true_SumHEM[[simID]] <- con_mat(b1 = df_est$b_true, b2 = df_est$b_SumHEM)
        }
        # f1-score of b_est vs true
        {
          ls_f1_b_true_true[[simID]] <- f1_score(mat = ls_mat_b_true_true[[simID]])
          ls_f1_b_true_gwas[[simID]] <- f1_score(mat = ls_mat_b_true_gwas[[simID]])
          ls_f1_b_true_infLDpred[[simID]] <- f1_score(mat = ls_mat_b_true_infLDpred[[simID]])
          ls_f1_b_true_autoLDpred[[simID]] <- f1_score(mat = ls_mat_b_true_autoLDpred[[simID]])
          ls_f1_b_true_SumRR[[simID]] <- f1_score(mat = ls_mat_b_true_SumRR[[simID]])
          ls_f1_b_true_SumHEM[[simID]] <- f1_score(mat = ls_mat_b_true_SumHEM[[simID]])
        }
        # dataframe
        {
          df <- data.frame(simID = simID, Nsim = Nsim, N_train = N_train, N_test = N_test, M = M, h2 = h2, pi = pi, Mcausal = Mcausal, h2input = h2input, ws = ws, 
                           macrof1_true = ls_f1_b_true_true[[simID]]$macro_f1, 
                           macrof1_gwas = ls_f1_b_true_gwas[[simID]]$macro_f1, 
                           macrof1_inf_LDpred = ls_f1_b_true_infLDpred[[simID]]$macro_f1, 
                           macrof1_auto_LDpred = ls_f1_b_true_autoLDpred[[simID]]$macro_f1, 
                           macrof1_SumRR = ls_f1_b_true_SumRR[[simID]]$macro_f1, 
                           macrof1_SumHEM = ls_f1_b_true_SumHEM[[simID]]$macro_f1)
        }
        df_eval_f1 <- rbind(df_eval_f1, df)
      }
      # save
      {
        if (h2input==h2 & ws==500) {
          saveRDS(df_eval_f1, file = paste0(workpath,"/df_eval_f1.rds"))
          saveRDS(list(ls_mat_b_true_true = ls_mat_b_true_true, ls_mat_b_true_gwas = ls_mat_b_true_gwas, 
                       ls_mat_b_true_infLDpred = ls_mat_b_true_infLDpred, ls_mat_b_true_autoLDpred = ls_mat_b_true_autoLDpred, 
                       ls_mat_b_true_SumRR = ls_mat_b_true_SumRR, ls_mat_b_true_SumHEM = ls_mat_b_true_SumHEM, 
                       ls_f1_b_true_true = ls_f1_b_true_true, ls_f1_b_true_gwas = ls_f1_b_true_gwas, 
                       ls_f1_b_true_infLDpred = ls_f1_b_true_infLDpred, ls_f1_b_true_autoLDpred = ls_f1_b_true_autoLDpred, 
                       ls_f1_b_true_SumRR = ls_f1_b_true_SumRR, ls_f1_b_true_SumHEM = ls_f1_b_true_SumHEM), 
                  file = paste0(workpath,"/ls_eval_conmat_f1.rds"))
        } else if (h2input!=h2 & ws==500) {
          saveRDS(df_eval_f1, file = paste0(workpath,"/df_eval_f1_h2input",h2input,".rds"))
          saveRDS(list(ls_mat_b_true_true = ls_mat_b_true_true, ls_mat_b_true_gwas = ls_mat_b_true_gwas, 
                       ls_mat_b_true_infLDpred = ls_mat_b_true_infLDpred, ls_mat_b_true_autoLDpred = ls_mat_b_true_autoLDpred, 
                       ls_mat_b_true_SumRR = ls_mat_b_true_SumRR, ls_mat_b_true_SumHEM = ls_mat_b_true_SumHEM, 
                       ls_f1_b_true_true = ls_f1_b_true_true, ls_f1_b_true_gwas = ls_f1_b_true_gwas, 
                       ls_f1_b_true_infLDpred = ls_f1_b_true_infLDpred, ls_f1_b_true_autoLDpred = ls_f1_b_true_autoLDpred, 
                       ls_f1_b_true_SumRR = ls_f1_b_true_SumRR, ls_f1_b_true_SumHEM = ls_f1_b_true_SumHEM), 
                  file = paste0(workpath,"/ls_eval_conmat_f1_h2input",h2input,".rds"))
        } else if (h2input==h2 & ws!=500) {
          saveRDS(df_eval_f1, file = paste0(workpath,"/df_eval_f1_ws",ws,".rds"))
          saveRDS(list(ls_mat_b_true_true = ls_mat_b_true_true, ls_mat_b_true_gwas = ls_mat_b_true_gwas, 
                       ls_mat_b_true_infLDpred = ls_mat_b_true_infLDpred, ls_mat_b_true_autoLDpred = ls_mat_b_true_autoLDpred, 
                       ls_mat_b_true_SumRR = ls_mat_b_true_SumRR, ls_mat_b_true_SumHEM = ls_mat_b_true_SumHEM, 
                       ls_f1_b_true_true = ls_f1_b_true_true, ls_f1_b_true_gwas = ls_f1_b_true_gwas, 
                       ls_f1_b_true_infLDpred = ls_f1_b_true_infLDpred, ls_f1_b_true_autoLDpred = ls_f1_b_true_autoLDpred, 
                       ls_f1_b_true_SumRR = ls_f1_b_true_SumRR, ls_f1_b_true_SumHEM = ls_f1_b_true_SumHEM), 
                  file = paste0(workpath,"/ls_eval_conmat_f1_ws",ws,".rds"))
        }
      }
      # }
    # }
  }
}

##### evaluation: beta distribution
ind_train = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_train.rds"))
ind_test = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_test.rds"))
nSNP = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_map/nSNP.rds"))
N_train = length(ind_train)
N_test = length(ind_test)
M = sum(nSNP)
require(dplyr)
h2 = 0.5
pi = 0.7
h2input = h2
ws = 500
Nsim = 20
n_bins = 32
for (h2 in c(0.1,0.5,0.7)) {
  for (pi in c(0.1,0.3,0.5,0.7)) {
    # for (h2input in c(0.3,0.4,0.5,0.6,0.7)) {
      # for (ws in c(250,500,1000)) {
      workpath = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/continuous/","h2",h2,"_pi",pi)
      Mcausal = round(M*pi)
      h2input = h2
      df_eval_corbins = df_eval_medbins = df_eval_meanbins <- NULL
      for (simID in 1:Nsim) {
        print(paste0("===== evaluation: beta distribution | h2=",h2," pi=",pi," h2input=",h2input," ws=",ws," simID=",simID,"/",Nsim," ====="))
        # load the data
        {
          if (h2input==h2 & ws==500) {
            df_est <- readRDS(paste0(workpath,"/res_est/df_est_simID",simID,".rds"))
          } else if (h2input!=h2 & ws==500) {
            df_est <- readRDS(paste0(workpath,"/res_est_h2input",h2input,"/df_est_simID",simID,".rds"))
          } else if (h2input==h2 & ws!=500) {
            df_est <- readRDS(paste0(workpath,"/res_est_ws",ws,"/df_est_simID",simID,".rds"))
          }
        }
        bins <- seq(from = 0, to = 1, by = 1/n_bins)
        qt <- quantile(df_est$b_true, probs = bins)
        # beta in bins
        df_betadist <- NULL
        {
          df_betadist <- df_est %>% mutate(simID = simID, Nsim = Nsim, N_train = N_train, N_test = N_test, M = M, h2 = h2, pi = pi, Mcausal = Mcausal, h2input = h2input, ws = ws, 
                                           n_bins = n_bins, bins = NA)
          for (k in 1:n_bins) {
            if (k==1) {
              bin_st = qt[k]-1e-10
            } else {
              bin_st = qt[k]
            }
            bin_ed = qt[k+1]
            df_betadist$bins[df_est$b_true>bin_st & df_est$b_true<=bin_ed] <- paste0("bin_",k)
          }
          
        }
        # save beta
        {
          # if (h2input==h2 & ws==500) {
          #   if (!fs::dir_exists(paste0(workpath,"/res_betadist"))) fs::dir_create(paste0(workpath,"/res_betadist"))
          #   saveRDS(df_betadist,
          #           file = paste0(workpath,"/res_betadist","/df_betadist_simID",simID,".rds"))
          # } else if (h2input!=h2 & ws==500) {
          #   if (!fs::dir_exists(paste0(workpath,"/res_betadist_h2input",h2input))) fs::dir_create(paste0(workpath,"/res_betadist_h2input",h2input))
          #   saveRDS(df_betadist,
          #           file = paste0(workpath,"/res_betadist_h2input",h2input,"/df_betadist_simID",simID,".rds"))
          # } else if (h2input==h2 & ws!=500) {
          #   if (!fs::dir_exists(paste0(workpath,"/res_betadist_ws",ws))) fs::dir_create(paste0(workpath,"/res_betadist_ws",ws))
          #   saveRDS(df_betadist,
          #           file = paste0(workpath,"/res_betadist_ws",ws,"/df_betadist_simID",simID,".rds"))
          # }
        }
        # statistic in bins
        df_corbins = df_medbins = df_meanbins <- NULL
        {
          numbins_true <- NULL
          meanbins_true = meanbins_gwas = meanbins_inf_LDpred = meanbins_auto_LDpred = meanbins_SumRR = meanbins_SumHEM <- numeric(n_bins)
          for (k in 1:n_bins) {
            if (k==1) {
              bin_st = qt[k]-1e-10
            } else {
              bin_st = qt[k]
            }
            bin_ed = qt[k+1]
            numbins_true[k] <- sum(df_est$b_true>bin_st & df_est$b_true<=bin_ed)
            
            meanbins_true[k] <- mean(df_est$b_true[df_est$b_true>bin_st & df_est$b_true<=bin_ed])
            meanbins_gwas[k] <- mean(df_est$sbeta[df_est$b_true>bin_st & df_est$b_true<=bin_ed])
            meanbins_inf_LDpred[k] <- mean(df_est$b_inf_LDpred[df_est$b_true>bin_st & df_est$b_true<=bin_ed])
            meanbins_auto_LDpred[k] <- mean(df_est$b_auto_LDpred[df_est$b_true>bin_st & df_est$b_true<=bin_ed])
            meanbins_SumRR[k] <- mean(df_est$b_SumRR[df_est$b_true>bin_st & df_est$b_true<=bin_ed])
            meanbins_SumHEM[k] <- mean(df_est$b_SumHEM[df_est$b_true>bin_st & df_est$b_true<=bin_ed])
            
          }
          df_meanbins <- as.data.frame(matrix(c(numbins_true,meanbins_true,meanbins_gwas,meanbins_inf_LDpred,meanbins_auto_LDpred,meanbins_SumRR,meanbins_SumHEM), nrow = 7, byrow = TRUE))
          colnames(df_meanbins) <- paste0("bin_",1:n_bins)
          df_meanbins <- df_meanbins %>% mutate(simID = simID, Nsim = Nsim, N_train = N_train, N_test = N_test, M = M, h2 = h2, pi = pi, Mcausal = Mcausal, h2input = h2input, ws = ws, 
                                                statistic = "meanbins", 
                                                method = c("nSNP in bins","True","GWAS","LDpred2-inf","LDpred2-auto","SumRR","SumHEM"))
          
        }
        df_eval_meanbins <- rbind(df_eval_meanbins, df_meanbins)
      }
      # save statistic
      {
        if (h2input==h2 & ws==500) {
          saveRDS(list(df_eval_meanbins = df_eval_meanbins), 
                  file = paste0(workpath,"/ls_eval_betadist.rds"))
        } else if (h2input!=h2 & ws==500) {
          saveRDS(list(df_eval_meanbins = df_eval_meanbins), 
                  file = paste0(workpath,"/ls_eval_betadist_h2input",h2input,".rds"))
        } else if (h2input==h2 & ws!=500) {
          saveRDS(list(df_eval_meanbins = df_eval_meanbins), 
                  file = paste0(workpath,"/ls_eval_betadist_ws",ws,".rds"))
        }
      }
      # }
    # }
  }
}

##### evaluation: se of h2_SumHEM
ind_train = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_train.rds"))
ind_test = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_test.rds"))
nSNP = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_map/nSNP.rds"))
N_train = length(ind_train)
N_test = length(ind_test)
M = sum(nSNP)
require(dplyr)
ws = 500
Nsim = 20
# h2 = 0.5
# pi = 0.1
for (h2 in c(0.1,0.5,0.7)) {
  for (pi in c(0.1,0.3,0.5,0.7)) {
    workpath = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/continuous/","h2",h2,"_pi",pi)
    Mcausal = round(M*pi)
    h2input = h2
    df_eval_h2se <- NULL
    for (simID in 1:Nsim) {
      print(paste0("===== evaluation: se of h2_SumHEM | h2=",h2," pi=",pi," simID=",simID,"/",Nsim," ====="))
      # load the data
      {
        if (h2input==h2) {
          df_est <- readRDS(paste0(workpath,"/res_est/df_est_simID",simID,".rds"))
        } else {
          df_est <- readRDS(paste0(workpath,"/res_est_h2input",h2input,"/df_est_simID",simID,".rds"))
        }
        df_h2var <- df_est %>% select(chr,winID,h2var_SumHEM) %>% distinct(chr,winID, .keep_all = TRUE)
      }
      # dataframe
      {
        df <- data.frame(simID = simID, Nsim = Nsim, N_train = N_train, N_test = N_test, M = M, h2 = h2, pi = pi, Mcausal = Mcausal, h2input = h2input, ws = ws, 
                         h2se_SumHEM = sqrt(sum(df_h2var$h2var_SumHEM)))
      }
      df_eval_h2se <- rbind(df_eval_h2se, df)
    }
    # save
    {
      if (h2input==h2) {
        saveRDS(df_eval_h2se, file = paste0(workpath,"/df_eval_h2se.rds"))
      } else {
        saveRDS(df_eval_h2se, file = paste0(workpath,"/df_eval_h2se_h2input",h2input,".rds"))
      }
    }
  }
}

##### evaluation: h2 (LDpred2-auto without h2 updating)
ind_train = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_train.rds"))
ind_test = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/ind_test.rds"))
nSNP = readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_map/nSNP.rds"))
N_train = length(ind_train)
N_test = length(ind_test)
M = sum(nSNP)
require(dplyr)
h2 = 0.5
pi = 0.7
h2input = h2
ws = 500
Nsim = 20
for (h2 in c(0.5)) {
  for (pi in c(0.1,0.3,0.5,0.7)) {
    workpath = paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/continuous/","h2",h2,"_pi",pi)
    Mcausal = round(M*pi)
    # h2input = h2
    df_eval_h2 <- NULL
    for (simID in 1:Nsim) {
      print(paste0("===== evaluation: h2 (LDpred-my) | h2=",h2," pi=",pi," h2input=",h2input," ws=",ws," simID=",simID,"/",Nsim," ====="))
      # load the data
      {
        df_est <- readRDS(paste0(workpath,"/res_est_myLDpred/df_est_simID",simID,".rds"))
      }
      # dataframe
      {
        df <- data.frame(simID = simID, Nsim = Nsim, N_train = N_train, N_test = N_test, M = M, h2 = h2, pi = pi, Mcausal = Mcausal, h2input = h2input, ws = ws, 
                         h2_true = as.numeric(df_est$b_true %*% df_est$Rb_true),
                         h2_gwas = as.numeric(df_est$sbeta %*% df_est$Rb_gwas),
                         h2_inf_LDpred = as.numeric(df_est$b_inf_LDpred %*% df_est$Rb_inf_LDpred), 
                         h2_auto_LDpred = as.numeric(df_est$b_auto_LDpred %*% df_est$Rb_auto_LDpred), 
                         h2_my_LDpred = as.numeric(df_est$b_my_LDpred %*% df_est$Rb_my_LDpred), 
                         h2_SumRR = as.numeric(df_est$b_SumRR %*% df_est$Rb_SumRR), 
                         h2_SumHEM = as.numeric(df_est$b_SumHEM %*% df_est$Rb_SumHEM))
      }
      df_eval_h2 <- rbind(df_eval_h2, df)
    }
    # save
    {
      saveRDS(df_eval_h2, file = paste0(workpath,"/df_eval_h2_myLDpred.rds"))
    }
  }
}







