###############################################################################
# Project: SumHEM
# Task: UKB validation | estimate SNP effects of phenotype
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
phenoID <- gsub(x = args[grep(x = args, pattern = "phenoID=")], pattern = "phenoID=", replacement = "")
file_est <- gsub(x = args[grep(x = args, pattern = "file_est=")], pattern = "file_est=", replacement = "")
h2 <- gsub(x = args[grep(x = args, pattern = "h2=")], pattern = "h2=", replacement = "")
ws <- gsub(x = args[grep(x = args, pattern = "ws=")], pattern = "ws=", replacement = "")
NCORES <- gsub(x = args[grep(x = args, pattern = "NCORES=")], pattern = "NCORES=", replacement = "")
path_out <- gsub(x = args[grep(x = args, pattern = "path_out=")], pattern = "path_out=", replacement = "")

h2 <- as.numeric(h2)
ws <- as.numeric(ws)
NCORES <- as.numeric(NCORES)

require(dplyr)
require(Matrix)
require(doParallel)
tryCatch(expr = {
  st = Sys.time()
  # GWAS
  cat(paste0("(",phenoID,") Loading GWAS data ...\n"))
  df_gwas <- readRDS(file_est) %>%
    select(rsid, chr, pos, beta, se, n_complete_samples) %>%
    rename(SNP = rsid, n_eff = n_complete_samples, beta_se = se) %>%
    mutate(na_SNP = is.na(beta), 
           beta = ifelse(is.na(beta),0,beta), 
           beta_se = ifelse(is.na(beta_se),1,beta_se), 
           n_eff = ifelse(is.na(n_eff),median(n_eff,na.rm = TRUE),n_eff), 
           scale = sqrt(beta^2+n_eff*beta_se^2), 
           sbeta = beta/scale)
  # LD info
  cat(paste0("(",phenoID,") Loading LD info ...\n"))
  df_map <- readRDS("/opt/working/projects/prj_035_LDpredRR/1000G/LD_EUR/map.rds") %>% select(rsid, chr, pos) %>% rename(SNP = rsid)
  # function
  estimation_ws <- function(j, ws, df_gwas_chr, df_map_chr, R_chr, h2, M) {
    print(paste0("estimation winID ", j, "/", ceiling(nrow(df_gwas_chr)/ws)))
    st_BP <- ws*(j-1) + 1
    ed_BP <- ifelse(ws*j>=(nrow(df_gwas_chr)-round(ws/2)), nrow(df_gwas_chr), ws*j)
    if (j!=1 & (ed_BP-st_BP+1) < (ws-round(ws/2))) return(NULL)
    df_gwas_ws <- df_gwas_chr[st_BP:ed_BP, ]
    R_ws <- R_chr[df_map_chr$mapID[st_BP:ed_BP], df_map_chr$mapID[st_BP:ed_BP]]
    R_SFBM <- bigsparser::as_SFBM(R_ws)
    res_ws <- list(h2 = h2, winSize = ncol(R_ws), winID = j)
    print("LDpred-inf")
    {
      res_ws$b_inf_LDpred <- bigsnpr::snp_ldpred2_inf(corr = R_SFBM, df_beta = df_gwas_ws, h2 = h2*ncol(R_SFBM)/M) / df_gwas_ws$scale
      res_ws$h2_inf_LDpred <- as.numeric(res_ws$b_inf_LDpred %*% R_ws %*% res_ws$b_inf_LDpred)
      res_ws$Rb_inf_LDpred <- as.numeric(R_ws %*% res_ws$b_inf_LDpred)
    }
    print("LDpred-auto")
    {
      model <- bigsnpr::snp_ldpred2_auto(corr = R_SFBM, df_beta = df_gwas_ws, 
                                         h2_init = h2*ncol(R_SFBM)/M, ncores = 1, 
                                         vec_p_init = bigsnpr::seq_log(1e-4,0.2,length.out=30), 
                                         allow_jump_sign = FALSE, shrink_corr = 0.95)
      range <- sapply(model, function(auto) diff(range(auto$corr_est)))
      keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))
      res_ws$b_auto_LDpred <- rowMeans(sapply(model[keep], function(auto) auto$beta_est)) / df_gwas_ws$scale
      res_ws$pp_auto_LDpred <- rowMeans(sapply(model[keep], function(auto) auto$postp_est))
      res_ws$pi_auto_LDpred <- mean(sapply(model[keep], function(auto) auto$p_est))
      res_ws$h2_auto_LDpred <- as.numeric(res_ws$b_auto_LDpred %*% R_ws %*% res_ws$b_auto_LDpred)
      res_ws$Rb_auto_LDpred <- as.numeric(R_ws %*% res_ws$b_auto_LDpred)
      rm(model)
    }
    print("SumRR & SumHEM")
    {
      # SumRR
      res_ws$sig2e <- 1-h2*ncol(R_SFBM)/M
      res_ws$sig2b_SumRR <- h2/M
      res_ws$lam_SumRR <- res_ws$sig2e / res_ws$sig2b_SumRR
      res_ws$b_SumRR <- bigsparser::sp_solve_sym(A = R_SFBM, b = df_gwas_ws$sbeta, add_to_diag = res_ws$lam_SumRR/df_gwas_ws$n_eff)
      res_ws$h2_SumRR <- as.numeric(res_ws$b_SumRR %*% R_ws %*% res_ws$b_SumRR)
      res_ws$Rb_SumRR <- as.numeric(R_ws %*% res_ws$b_SumRR)
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
      rm(Q)
      
    }
    return(list(res_ws = res_ws))
  }
  # estimates
  df_est <- NULL
  for (CHR in 1:22) {
    st_chr = Sys.time()
    cat(paste0("(",phenoID,") Estimating with SumHEM CHR ",CHR,"/22 ...\t"))
    # data
    df_gwas_chr <- df_gwas %>% filter(chr == CHR)
    df_map_chr <- df_map %>% filter(chr == CHR) %>% mutate(mapID = rank(pos)) %>% right_join(., df_gwas_chr, by = c("SNP","chr","pos"))
    R_chr <- readRDS(paste0("/opt/working/projects/prj_035_LDpredRR/1000G/LD_EUR/LD_chr", CHR, ".rds"))
    R_chr@x[is.na(R_chr@x)] <- 1e-10
    # computation
    cl <- makeCluster(NCORES)
    registerDoParallel(cl)
    res_chr <- foreach::foreach(
      j = 1:ceiling(nrow(df_gwas_chr)/ws),
      .combine = c, 
      .packages = c("bigsnpr", "Matrix", "bigsparser", "bigstatsr", "dplyr")
    ) %dopar% estimation_ws(j = j, ws = ws, df_gwas_chr = df_gwas_chr, df_map_chr = df_map_chr, R_chr = R_chr, 
                            h2 = h2, M = nrow(df_gwas))
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
    cat(paste0("Time ",difftime(ed_chr,st_chr,units="mins")," min \n"))
  }
  # save
  saveRDS(df_est, file = paste0(path_out,"/df_est.",phenoID,".rds"))
  rm(df_gwas, df_est)
  ed = Sys.time()
  cat(paste0("(",phenoID,") Estimation is Done !!! Time ",difftime(ed,st,units="mins")," min \n"))
}, error = function(e){
  cat(e)
  ed = Sys.time()
  stop(paste0("(",phenoID,") Error !!! Time ",difftime(ed,st,units="mins")," min \n"))
})
























