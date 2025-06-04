###############################################################################
# Project: SumHEM
# Task: simulation | simulate SNP effects and phenotype
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
h2 <- gsub(x = args[grep(x = args, pattern = "h2=")], pattern = "h2=", replacement = "")
pi <- gsub(x = args[grep(x = args, pattern = "pi=")], pattern = "pi=", replacement = "")
Mcausal <- gsub(x = args[grep(x = args, pattern = "Mcausal=")], pattern = "Mcausal=", replacement = "")
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
if (length(pi)!=0) {
  pi <- as.numeric(pi)
  Mcausal <- round(M*pi)
} else if (length(Mcausal)!=0) {
  Mcausal <- as.numeric(Mcausal)
  pi <- Mcausal/M
}
ws <- as.numeric(ws)
simID <- as.numeric(simID)

##### create direction
{
  if (!fs::dir_exists(workpath)) fs::dir_create(workpath)
  if (!fs::dir_exists(paste0(workpath,"/res_true"))) fs::dir_create(paste0(workpath,"/res_true"))
  if (!fs::dir_exists(paste0(workpath,"/res_gwas"))) fs::dir_create(paste0(workpath,"/res_gwas"))
  if (!fs::dir_exists(paste0(workpath,"/res_h2"))) fs::dir_create(paste0(workpath,"/res_h2"))
}

##### simulation
st_ID = Sys.time()
{
  cat(paste0("========== simulation h2=",h2," pi=",round(pi,4)," Mcausal=",Mcausal," simID=",simID," ==========","\n"))
  # true SNP effects
  b_true <- numeric(M)
  ind_causal <- sort(sample(1:M, size = Mcausal))
  b_true[ind_causal] <- rnorm(n = Mcausal, mean = 0, sd = sqrt(h2/Mcausal))
  df_beta <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_map/df_map.rds")) %>%
    mutate(b_true = b_true)
  # genetic components
  y_gen <- numeric(N_train+N_test)
  for (CHR in 1:22) {
    cat(paste0("simulation CHR ", CHR, "/22","\t"))
    st_chr = Sys.time()
    obj.bigSNP <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_imputed/chr_",CHR,".rds"))
    obj.bigSNP$genotypes$backingfile <- paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_imputed/chr_",CHR,".bk")
    obj.bigSNP$genotypes$code256 <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_imputed_codes/chr_",CHR,".rds"))$code256
    sG_train <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_scaled/sG_train_chr", CHR, ".rds"))
    sG_test <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_scaled/sG_test_chr", CHR, ".rds"))
    df_beta_chr <- df_beta %>% filter(chr == CHR)
    y_gen[ind_train] <- y_gen[ind_train] + 
      bigstatsr::big_prodVec(X = obj.bigSNP$genotypes, ind.row = ind_train, y.col = df_beta_chr$b_true, 
                             center = sG_train$center, scale = sG_train$scale, ncores = 1)
    y_gen[ind_test] <- y_gen[ind_test] + 
      bigstatsr::big_prodVec(X = obj.bigSNP$genotypes, ind.row = ind_test, y.col = df_beta_chr$b_true, 
                             center = sG_test$center, scale = sG_test$scale, ncores = 1)
    rm(obj.bigSNP,sG_train,sG_test,df_beta_chr)
    ed_chr = Sys.time()
    cat(paste0("Done !!! Time ",difftime(ed_chr,st_chr,units="mins")," min \n"))
  }
  # phenotype
  y_env <- rnorm(n = length(y_gen), mean = 0, sd = sqrt(1-h2))
  y_sim <- y_gen + y_env
  # save
  saveRDS(ind_causal, file = paste0(workpath,"/res_true/ind_causal_simID",simID,".rds"))
  saveRDS(b_true, file = paste0(workpath,"/res_true/b_true_simID",simID,".rds"))
  saveRDS(y_sim, file = paste0(workpath,"/res_true/y_sim_simID",simID,".rds"))
  saveRDS(y_gen, file = paste0(workpath,"/res_true/y_gen_simID",simID,".rds"))
  saveRDS(y_env, file = paste0(workpath,"/res_true/y_env_simID",simID,".rds"))
  
  rm(df_beta,ind_causal,b_true,y_sim,y_gen,y_env)
}
ed_ID = Sys.time()
cat(paste0("Simulation is done !!! Time ",difftime(ed_ID,st_ID,units="mins")," min \n"))

##### GWAS
st_ID = Sys.time()
{
  cat(paste0("========== GWAS h2=",h2," pi=",round(pi,4)," Mcausal=",Mcausal," simID=",simID," ==========","\n"))
  # load the data
  y_sim <- readRDS(paste0(workpath,"/res_true/y_sim_simID",simID,".rds"))
  # gwas
  b_gwas_train = b_gwas_se_train = b_gwas_test = b_gwas_se_test = NULL
  for (CHR in 1:22) {
    cat(paste0("GWAS CHR ", CHR, "/22","\t"))
    st_chr = Sys.time()
    obj.bigSNP <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_imputed/chr_",CHR,".rds"))
    obj.bigSNP$genotypes$backingfile <- paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_imputed/chr_",CHR,".bk")
    obj.bigSNP$genotypes$code256 <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_imputed_codes/chr_",CHR,".rds"))$code256
    # train
    df_train <- bigstatsr::big_univLinReg(X = obj.bigSNP$genotypes, y.train = y_sim[ind_train], ind.train = ind_train, ncores = 1)
    b_gwas_train <- c(b_gwas_train, ifelse(is.na(df_train$estim), 0, df_train$estim))
    b_gwas_se_train <- c(b_gwas_se_train, ifelse(is.na(df_train$std.err), 1, df_train$std.err))
    # test
    df_test <- bigstatsr::big_univLinReg(X = obj.bigSNP$genotypes, y.train = y_sim[ind_test], ind.train = ind_test, ncores = 1)
    b_gwas_test <- c(b_gwas_test, ifelse(is.na(df_test$estim), 0, df_test$estim))
    b_gwas_se_test <- c(b_gwas_se_test, ifelse(is.na(df_test$std.err), 1, df_test$std.err))
    
    rm(obj.bigSNP,df_train,df_test)
    ed_chr = Sys.time()
    cat(paste0("Done !!! Time ",difftime(ed_chr,st_chr,units="mins")," min \n"))
  }
  # save
  saveRDS(list(b_gwas = b_gwas_train, b_gwas_se = b_gwas_se_train), 
          file = paste0(workpath,"/res_gwas/ls_gwas_train_simID",simID,".rds"))
  saveRDS(list(b_gwas = b_gwas_test, b_gwas_se = b_gwas_se_test), 
          file = paste0(workpath,"/res_gwas/ls_gwas_test_simID",simID,".rds"))
  
  rm(y_sim,b_gwas_train,b_gwas_se_train,b_gwas_test,b_gwas_se_test)
}
ed_ID = Sys.time()
cat(paste0("GWAS is done !!! Time ",difftime(ed_ID,st_ID,units="mins")," min \n"))

##### LDSC heritability
st_ID = Sys.time()
{
  cat(paste0("========== heritability h2=",h2," pi=",round(pi,4)," Mcausal=",Mcausal," simID=",simID," ==========","\n"))
  # load the data
  LDscore <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/R_dsCM/df_LD.rds"))$LDscore
  res_gwas <- readRDS(paste0(workpath,"/res_gwas/ls_gwas_train_simID",simID,".rds"))
  df_map <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_map/df_map.rds"))
  # LDSC h2
  {
    df_gwas <- df_map %>%
      mutate(beta = res_gwas$b_gwas, beta_se = res_gwas$b_gwas_se, n_eff = N_train)
    res_ldsc <- with(df_gwas, bigsnpr::snp_ldsc(LDscore, length(LDscore), chi2 = (beta/beta_se)^2, sample_size = n_eff, blocks = NULL))
  }
  df_h2 <- data.frame(simID = simID, N = N_train, M = M, Mcausal = Mcausal, pi = pi, h2 = h2, h2_ldsc = res_ldsc[["h2"]])
  # save
  saveRDS(df_h2, file = paste0(workpath,"/res_h2/df_h2_simID",simID,".rds"))
  
  rm(LDscore,res_gwas,df_map,df_gwas,res_ldsc,df_h2)
}
ed_ID = Sys.time()
cat(paste0("Heritability computation is done !!! Time ",difftime(ed_ID,st_ID,units="mins")," min \n"))




























