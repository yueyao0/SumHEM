###############################################################################
# Project: SumHEM
# Task: simulation | estimate LDSC heritability and SNP effects
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
h2 <- gsub(x = args[grep(x = args, pattern = "h2=")], pattern = "h2=", replacement = "")
pi <- gsub(x = args[grep(x = args, pattern = "pi=")], pattern = "pi=", replacement = "")
Mcausal <- gsub(x = args[grep(x = args, pattern = "Mcausal=")], pattern = "Mcausal=", replacement = "")
piinput <- gsub(x = args[grep(x = args, pattern = "piinput=")], pattern = "piinput=", replacement = "")
h2input <- gsub(x = args[grep(x = args, pattern = "h2input=")], pattern = "h2input=", replacement = "")
ws <- gsub(x = args[grep(x = args, pattern = "ws=")], pattern = "ws=", replacement = "")
simID <- gsub(x = args[grep(x = args, pattern = "simID=")], pattern = "simID=", replacement = "")
workpath <- gsub(x = args[grep(x = args, pattern = "workpath=")], pattern = "workpath=", replacement = "")
NCORES <- gsub(x = args[grep(x = args, pattern = "NCORES=")], pattern = "NCORES=", replacement = "")
Mtopmiss <- gsub(x = args[grep(x = args, pattern = "Mtopmiss=")], pattern = "Mtopmiss=", replacement = "")

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
if (length(h2input)!=0) {
  h2input <- as.numeric(h2input)
} else {
  h2input <- h2
}
if (length(piinput)!=0) {
  piinput <- as.numeric(piinput)
} else {
  piinput <- pi
}
ws <- as.numeric(ws)
simID <- as.numeric(simID)
NCORES <- as.numeric(NCORES)
if (length(Mtopmiss)==0) {
  Mtopmiss <- 0
} else {
  Mtopmiss <- as.numeric(Mtopmiss)
}

##### create direction
{
  if (Mtopmiss==0 & h2input==h2 & ws==500) {
    if (!fs::dir_exists(paste0(workpath,"/res_h2"))) fs::dir_create(paste0(workpath,"/res_h2"))
    if (!fs::dir_exists(paste0(workpath,"/res_est"))) fs::dir_create(paste0(workpath,"/res_est"))
  } else if (Mtopmiss!=0 & h2input==h2 & ws==500) {
    if (!fs::dir_exists(paste0(workpath,"/res_h2_Mtopmiss",Mtopmiss))) fs::dir_create(paste0(workpath,"/res_h2_Mtopmiss",Mtopmiss))
    if (!fs::dir_exists(paste0(workpath,"/res_est_Mtopmiss",Mtopmiss))) fs::dir_create(paste0(workpath,"/res_est_Mtopmiss",Mtopmiss))
  } else if (Mtopmiss==0 & h2input!=h2 & ws==500) {
    if (!fs::dir_exists(paste0(workpath,"/res_h2_h2input",h2input))) fs::dir_create(paste0(workpath,"/res_h2_h2input",h2input))
    if (!fs::dir_exists(paste0(workpath,"/res_est_h2input",h2input))) fs::dir_create(paste0(workpath,"/res_est_h2input",h2input))
  } else if (Mtopmiss==0 & h2input==h2 & ws!=500) {
    if (!fs::dir_exists(paste0(workpath,"/res_h2_ws",ws))) fs::dir_create(paste0(workpath,"/res_h2_ws",ws))
    if (!fs::dir_exists(paste0(workpath,"/res_est_ws",ws))) fs::dir_create(paste0(workpath,"/res_est_ws",ws))
  }
}

##### load the data
{
  res_gwas <- readRDS(paste0(workpath,"/res_gwas/ls_gwas_train_simID",simID,".rds"))
  b_true <- readRDS(paste0(workpath,"/res_true/b_true_simID",simID,".rds"))
  df_map <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/G_map/df_map.rds"))
  LDscore <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/R_dsCM/df_LD.rds"))$LDscore
}

##### dataframe preparation & missing observations
{
  df_gwas <- df_map %>%
    mutate(beta = res_gwas$b_gwas, beta_se = res_gwas$b_gwas_se, n_eff = N_train, 
           scale = sqrt(beta^2+n_eff*beta_se^2), sbeta = beta/scale, 
           b_true = b_true)
  if (Mtopmiss!=0) {
    ind_miss <- order(abs(b_true), decreasing = TRUE)[1:Mtopmiss]
    M <- length(b_true)-length(ind_miss)
    df_gwas <- df_gwas[-ind_miss, ]
    LDscore <- LDscore[-ind_miss]
  }
}

##### heritability
st_ID = Sys.time()
{
  cat(paste0("========== heritability h2=",h2," pi=",round(pi,4)," Mcausal=",Mcausal," Mtopmiss=",Mtopmiss," h2input=",h2input," ws=",ws," simID=",simID," ==========","\n"))
  # LDSC h2
  res_ldsc <- with(df_gwas, bigsnpr::snp_ldsc(LDscore, length(LDscore), chi2 = (beta/beta_se)^2, sample_size = n_eff, blocks = NULL))
  df_h2 <- data.frame(simID = simID, N = N_train, M = M, Mcausal = Mcausal, pi = pi, Mtopmiss = Mtopmiss, h2 = h2, h2_ldsc = res_ldsc[["h2"]])
  # save
  if (Mtopmiss==0 & h2input==h2 & ws==500) {
    saveRDS(df_h2, file = paste0(workpath,"/res_h2/df_h2_simID",simID,".rds"))
  } else if (Mtopmiss!=0 & h2input==h2 & ws==500) {
    saveRDS(df_h2, file = paste0(workpath,"/res_h2_Mtopmiss",Mtopmiss,"/df_h2_simID",simID,".rds"))
  } else if (Mtopmiss==0 & h2input!=h2 & ws==500) {
    saveRDS(df_h2, file = paste0(workpath,"/res_h2_h2input",h2input,"/df_h2_simID",simID,".rds"))
  } else if (Mtopmiss==0 & h2input==h2 & ws!=500) {
    saveRDS(df_h2, file = paste0(workpath,"/res_h2_ws",ws,"/df_h2_simID",simID,".rds"))
  }
  rm(res_ldsc, df_h2)
}
ed_ID = Sys.time()
cat(paste0("Heritability computation is done !!! Time ",difftime(ed_ID,st_ID,units="mins")," min \n"))

##### estimation
estimation_ws <- function(j,ws,df_gwas_chr,df_map_chr,R_chr,pi_est,h2_est,M) {
  print(paste0("estimation winID ", j, "/", ceiling(nrow(df_gwas_chr)/ws)))
  st_BP <- ws*(j-1) + 1
  ed_BP <- ifelse(ws*j>=(nrow(df_gwas_chr)-round(ws/2)), nrow(df_gwas_chr), ws*j)
  if (j!=1 & (ed_BP-st_BP+1) < (ws-round(ws/2))) return(NULL)
  df_gwas_ws <- df_gwas_chr[st_BP:ed_BP, ]
  R_ws <- R_chr[df_map_chr$mapID[st_BP:ed_BP], df_map_chr$mapID[st_BP:ed_BP]]
  R_SFBM <- bigsparser::as_SFBM(R_ws)
  res_ws <- list(M = M, h2_est = h2_est, h2_est_ws = h2_est*nrow(df_gwas_ws)/M, 
                 winSize = ncol(R_ws), winID = j)
  print("LDpred-inf")
  st_ws = Sys.time()
  {
    res_ws$b_inf_LDpred <- bigsnpr::snp_ldpred2_inf(corr = R_SFBM, df_beta = df_gwas_ws, h2 = h2_est*nrow(df_gwas_ws)/M) / df_gwas_ws$scale
    res_ws$h2_inf_LDpred <- as.numeric(res_ws$b_inf_LDpred %*% R_ws %*% res_ws$b_inf_LDpred)
    res_ws$Rb_inf_LDpred <- as.numeric(R_ws %*% res_ws$b_inf_LDpred)
  }
  ed_ws = Sys.time()
  res_ws$time_inf_LDpred <- as.numeric(difftime(ed_ws,st_ws,units="mins"))
  
  print("LDpred-auto")
  st_ws = Sys.time()
  {
    model <- bigsnpr::snp_ldpred2_auto(corr = R_SFBM, df_beta = df_gwas_ws, vec_p_init = pi_est, 
                                       h2_init = h2_est*nrow(df_gwas_ws)/M, ncores = 1)
    res_ws$b_auto_LDpred <- model[[1]]$beta_est / df_gwas_ws$scale
    res_ws$pp_auto_LDpred <- model[[1]]$postp_est
    res_ws$pi_auto_LDpred <- model[[1]]$p_est
    res_ws$h2_auto_LDpred <- as.numeric(res_ws$b_auto_LDpred %*% R_ws %*% res_ws$b_auto_LDpred)
    res_ws$Rb_auto_LDpred <- as.numeric(R_ws %*% res_ws$b_auto_LDpred)
    rm(model)
  }
  ed_ws = Sys.time()
  res_ws$time_auto_LDpred <- as.numeric(difftime(ed_ws,st_ws,units="mins"))
  
  print("SumRR & SumHEM")
  {
    st_ws = Sys.time()
    # SumRR
    res_ws$sig2e <- 1-h2_est*nrow(df_gwas_ws)/M
    res_ws$sig2b_SumRR <- h2_est/M
    res_ws$lam_SumRR <- res_ws$sig2e / res_ws$sig2b_SumRR
    res_ws$b_SumRR <- bigsparser::sp_solve_sym(A = R_SFBM, b = df_gwas_ws$sbeta, add_to_diag = res_ws$lam_SumRR/df_gwas_ws$n_eff)
    res_ws$h2_SumRR <- as.numeric(res_ws$b_SumRR %*% R_ws %*% res_ws$b_SumRR)
    res_ws$Rb_SumRR <- as.numeric(R_ws %*% res_ws$b_SumRR)
    ed_ws = Sys.time()
    res_ws$time_SumRR <- as.numeric(difftime(ed_ws,st_ws,units="mins"))
    # SumHEM
    Q <- R_ws*median(df_gwas_ws$n_eff) + Matrix::.symDiagonal(ncol(R_ws), x = res_ws$lam_SumRR)
    res_ws$hv_SumHEM <- res_ws$lam_SumRR * diag(solve(Q))
    res_ws$hv_SumHEM <- ifelse(res_ws$hv_SumHEM>(1-1e-10), 1-1e-10, res_ws$hv_SumHEM)
    res_ws$sig2b_SumHEM <- res_ws$b_SumRR^2 / (1-res_ws$hv_SumHEM)
    res_ws$sig2b_SumHEM <- ifelse(res_ws$sig2b_SumHEM<1e-10, 1e-10, res_ws$sig2b_SumHEM)
    res_ws$lam_SumHEM <- res_ws$sig2e / res_ws$sig2b_SumHEM
    res_ws$b_SumHEM <- bigsparser::sp_solve_sym(A = R_SFBM, b = df_gwas_ws$sbeta, add_to_diag = res_ws$lam_SumHEM/df_gwas_ws$n_eff)
    res_ws$h2_SumHEM <- as.numeric(res_ws$b_SumHEM %*% R_ws %*% res_ws$b_SumHEM)
    res_ws$Rb_SumHEM <- as.numeric(R_ws %*% res_ws$b_SumHEM)
    ed_ws = Sys.time()
    res_ws$time_SumHEM <- as.numeric(difftime(ed_ws,st_ws,units="mins"))
    # variance of inferred h2
    H <- res_ws$lam_SumRR * solve(Q)
    S_2 <- res_ws$b_SumRR %*% t(res_ws$b_SumRR) / (1-H)
    S_2[S_2<1e-10 & S_2>0] <- 1e-10
    S_2[S_2>-1e-10 & S_2<0] <- -1e-10
    P_2 <- R_ws %*% S_2
    M_2 <- as.numeric(res_ws$b_SumHEM%*%R_ws%*%S_2%*%R_ws%*%res_ws$b_SumHEM)
    res_ws$h2var_SumHEM <- 2*sum(diag(P_2 %*% P_2)) + 4*M_2
    rm(Q,H,S_2,P_2,M_2)
  }
  
  print("True")
  st_ws = Sys.time()
  {
    res_ws$h2_true <- as.numeric(df_gwas_ws$b_true %*% R_ws %*% df_gwas_ws$b_true)
    res_ws$Rb_true <- as.numeric(R_ws %*% df_gwas_ws$b_true)
  }
  ed_ws = Sys.time()
  res_ws$time_true <- as.numeric(difftime(ed_ws,st_ws,units="mins"))
  
  print("GWAS")
  st_ws = Sys.time()
  {
    res_ws$h2_gwas <- as.numeric(df_gwas_ws$sbeta %*% R_ws %*% df_gwas_ws$sbeta)
    res_ws$Rb_gwas <- as.numeric(R_ws %*% df_gwas_ws$sbeta)
  }
  ed_ws = Sys.time()
  res_ws$time_gwas <- as.numeric(difftime(ed_ws,st_ws,units="mins"))
  
  return(list(res_ws = res_ws))
}
st_ID = Sys.time()
{
  cat(paste0("========== estimation h2=",h2," pi=",round(pi,4)," Mcausal=",Mcausal," Mtopmiss=",Mtopmiss," h2input=",h2input," ws=",ws," simID=",simID," ==========","\n"))
  # estimate
  df_est <- NULL
  for (CHR in 1:22) {
    st_chr = Sys.time()
    cat(paste0("estimation CHR ", CHR, "/22","\t"))
    # data
    df_gwas_chr <- df_gwas %>% filter(chr == CHR)
    df_map_chr <- df_map %>% filter(chr == CHR) %>% mutate(mapID = rank(pos)) %>% right_join(., df_gwas_chr, by = c("SNP","chr","pos","A1","A2"))
    R_chr <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/UKBB/simulation/R_dsCM/R_chr", CHR, ".rds"))
    R_chr@x[is.na(R_chr@x)] <- 1e-10
    # computation
    cl <- makeCluster(NCORES)
    registerDoParallel(cl)
    res_chr <- foreach::foreach(
      j = 1:ceiling(nrow(df_gwas_chr)/ws),
      .combine = c, 
      .packages = c("bigsnpr", "Matrix", "bigsparser", "bigstatsr", "dplyr")
    ) %dopar% estimation_ws(j = j, ws = ws, df_gwas_chr = df_gwas_chr, df_map_chr = df_map_chr, R_chr = R_chr, 
                            pi_est = piinput, h2_est = h2input, M = nrow(df_gwas))
    stopCluster(cl)
    # summary
    df_est_chr <- NULL
    for (j in 1:length(res_chr)) {
      df_est_chr <- rbind(df_est_chr, as.data.frame(res_chr[[j]]))
    }
    df_est_chr <- cbind(df_gwas_chr %>% rename(b_gwas = beta, b_gwas_se = beta_se, N = n_eff), df_est_chr)
    df_est <- rbind(df_est, df_est_chr)
    rm(df_gwas_chr, df_map_chr, R_chr, res_chr, df_est_chr)
    ed_chr = Sys.time()
    cat(paste0("Done !!! Time ",difftime(ed_chr,st_chr,units="mins")," min \n"))
  }
  # save
  if (Mtopmiss==0 & h2input==h2 & ws==500) {
    saveRDS(df_est, file = paste0(workpath,"/res_est/df_est_simID",simID,".rds"))
  } else if (Mtopmiss!=0 & h2input==h2 & ws==500) {
    saveRDS(df_est, file = paste0(workpath,"/res_est_Mtopmiss",Mtopmiss,"/df_est_simID",simID,".rds"))
  } else if (Mtopmiss==0 & h2input!=h2 & ws==500) {
    saveRDS(df_est, file = paste0(workpath,"/res_est_h2input",h2input,"/df_est_simID",simID,".rds"))
  } else if (Mtopmiss==0 & h2input==h2 & ws!=500) {
    saveRDS(df_est, file = paste0(workpath,"/res_est_ws",ws,"/df_est_simID",simID,".rds"))
  }
  rm(df_est)
}
ed_ID = Sys.time()
cat(paste0("Estimation is done !!! Time ",difftime(ed_ID,st_ID,units="mins")," min \n"))



























