#' Summary-level heteroscedastic effects model (SumHEM)
#' 
#' The function returns the estimates of SNP effects of one trait based on GWAS summary statistics. 
#' 
#' @param df_gwas A data frame including GWAS summary statistics of genetic variants for a trait. 
#' The input data frame should include following columns: 
#' rsid, SNP ID; chr, chromosome where the SNP is; pos, position of the SNP; 
#' beta, estimate of marginal effect in GWAS; se, standard error of the estimates of marginal effects in GWAS; N, sample size; 
#' @param df_map A data frame including SNP information of whole-genome linkage disequilibrium (LD) reference panel. 
#' The information should be totally matched with the SNP including in LD files and the input data frame should include following columns: 
#' rsid, SNP ID; chr, chromosome where the SNP is; pos, position of the SNP;
#' @param LD_path Path to the directory where LD information is stored. 
#' The LD matrix should be separated by chromosomes and stored in 'dsCMatrix' class as a 'rds' object in each LD file. 
#' The filenames should be 'LD_chr'+chromosome+'.rds' 
#' @param h2input Heritability of the trait estimated by external software. 
#' The input numeric value should be between 0 and 1. 
#' If missing, the value will be computed in progress. 
#' @param ws The number of SNP in a window for estimation. 
#' @param NCORES The number of cores using for parallel computation. 
#' @note Users can download the LD correlation matrices for European ancestry population, authored by Florian Privé. 
#' The download link can be found at https://figshare.com/articles/dataset/European_LD_reference/13034123
#' This is the LD reference (correlations between pairs of genetic variants) for 1,054,330 HapMap3 variants based on 362,320 European individuals of the UK biobank. 
#' 
#' @return A data frame is returned with:
#' \itemize{
#' \item{rsid }{SNP ID.}
#' \item{chr }{The chromosome where the SNP in.}
#' \item{pos }{The position of the SNP.}
#' \item{beta }{The estimate of marginal effect in GWAS.}
#' \item{se }{The standard error of beta.}
#' \item{N }{The sample size.}
#' \item{scale }{The scale used in marginal effects, scale=sqrt(beta^2+N*se^2).}
#' \item{sbeta }{The scaled marginal effect, sbeta=beta/scale.}
#' \item{h2input }{The pre-estimated heritability of trait.}
#' \item{winSize }{The number of SNP included in the window.}
#' \item{winID }{The window ID in each chromosome.}
#' \item{sig2e }{The (sigma_e)^2 used in the window.}
#' \item{sig2b_SumRR }{The (sigma_b)^2 of step1 used in the window.}
#' \item{lambda_SumRR }{The shrinkage parameter lambda of step1 used in the window.}
#' \item{b_SumRR }{SumRR estimates of genetic effects in the window.}
#' \item{Rb_SumRR }{LD %*% SumRR estimates in the window.}
#' \item{hv }{The hat values used in the window.}
#' \item{sig2b_SumHEM }{The (sigma_bj)^2 of step2 used in the window.}
#' \item{lambda_SumHEM }{The shrinkage parameter lambda of step2 used in the window.}
#' \item{b_SumHEM }{SumHEM estimates of genetic effects in the window.}
#' \item{Rb_SumHEM }{LD %*% SumHEM estimates in the window.}
#' \item{h2_SumHEM }{The inferred heritability of SumHEM in the window.}
#' \item{h2var_SumHEM }{The variance of SumHEM-inferred heritability.}
#' }
#' 
#' @author Yue Yao, Wenzhuo Lin, Xia Shen
#' 
#' @references 
#' ...
#' 
#' @seealso 
#' SumHEM tutorial on GitHub page: https://github.com/yueyao0/SumHEM
#' 
#' @examples 
#'\dontrun{

#' ## Load example GWAS summary statistics for standing height
#' df_gwas <- readRDS("SumHEM/data/gwas.height.rds")
#' ## Load SNP information of LD reference panel
#' df_map <- readRDS("SumHEM/data/map.rds")
#' 
#' ## Setting for parameters
#' ws <- 500 # the number of SNPs in a window for estimation
#' NCORES <- 10 # the number of cores using in parallel computation
#' LD_path <- "SumHEM/LD_EUR" # path of folder with LD files
#' h2input <- 0.4780 # heritability estimated by an external software (optional)
#' 
#' ## Estimate SNP effects via SumHEM
#' df_est <- SumHEM(
#'   df_gwas = df_gwas, df_map = df_map,
#'   h2input = h2input, # if missing, it will be computed during the process
#'   ws = ws, NCORES = NCORES, LD_path = LD_path
#' )
#' 
#' ## View the results
#' str(df_est)
#' 
#' ## Summaries the inferred heritability
#' df <- df_est %>% distinct(chr,winID, .keep_all = TRUE)
#' sum(df$h2_SumRR) # inferred heritability of SumRR
#' sum(df$h2_SumHEM) # inferred heritability of SumHEM
#' sqrt(sum(df$h2var_SumHEM)) # se of SumHEM-inferred heritability
#' 
#' }
#' @export
#' 


## SumHEM function
SumHEM <-
  function(df_gwas, df_map, LD_path, h2input="", ws=500, NCORES=1){
    # colnames(df_gwas) = c("rsid","chr","pos","beta","se","N")
    # colnames(df_map) = c("rsid","chr","pos","ld")
    ## dependencies
    require(dplyr)
    require(Matrix)
    time_start <- date()
    cat(paste0("Analysis started on ",time_start,"\n"))
    ## estimate LDSC heritability if missing "h2input"
    if (h2input=="") {
      st <- Sys.time()
      cat(paste0("Estimating LDSC heritability ...\t"))
      df <- df_gwas %>% inner_join(df_map, by = c("rsid","chr","pos"))
      LDscore <- df$ld
      res_ldsc <- with(df, 
                       bigsnpr::snp_ldsc(LDscore, length(LDscore), 
                                         chi2 = (beta/se)^2, sample_size = N, 
                                         blocks = NULL)
                       )
      h2input <- res_ldsc[["h2"]]
      ed <- Sys.time()
      cat(paste0("Time ",difftime(ed,st,units="mins")," min \n"))
      cat(paste0("LDSC estimates of heritability = ",round(h2input,8),"\n"))
    }
    if(is.na(h2input) | h2input<0 | h2input>1){
      error.message <- "The input of h2input has to be a number between 0 and 1. \n"
      stop(error.message)
    }
    ## SNP info in LD
    df_map <- df_map %>% select(rsid,chr,pos)
    ## GWAS
    df_gwas <- df_gwas %>%
      select(rsid,chr,pos,beta,se,N) %>%
      mutate(beta = ifelse(is.na(beta),0,beta), 
             se = ifelse(is.na(se),1,se), 
             N = ifelse(is.na(N),median(df_gwas$N,na.rm = TRUE),N), 
             scale = sqrt(beta^2+N*se^2), 
             sbeta = beta/scale) %>% 
      inner_join(df_map, by = c("rsid","chr","pos")) %>% 
      arrange(chr,pos)
    cat(paste0("The number of matched SNPs = ",nrow(df_gwas),"\n"))
    ## function: estimation in a window
    estimation_ws <- function(j, ws, df_gwas_chr, df_map_chr, R_chr, h2input, M){
      st_BP <- ws*(j-1) + 1
      ed_BP <- ifelse(ws*j>=(nrow(df_gwas_chr)-round(ws/2)), nrow(df_gwas_chr), ws*j)
      if (j!=1 & (ed_BP-st_BP+1) <= (ws-round(ws/2))) return(NULL)
      df_gwas_ws <- df_gwas_chr[st_BP:ed_BP, ]
      R_ws <- R_chr[df_map_chr$mapID[st_BP:ed_BP], df_map_chr$mapID[st_BP:ed_BP]]
      R_SFBM <- bigsparser::as_SFBM(R_ws)
      res_ws <- list(h2input = h2input, M = M, winSize = ncol(R_ws), winID = j)
      ## step1: SumRR estimates
      res_ws$sig2e <- 1-h2input*ncol(R_SFBM)/M
      res_ws$sig2b_SumRR <- h2input/M
      res_ws$lambda_SumRR <- res_ws$sig2e / res_ws$sig2b_SumRR
      res_ws$b_SumRR <- bigsparser::sp_solve_sym(A = R_SFBM, b = df_gwas_ws$sbeta, add_to_diag = res_ws$lambda_SumRR/df_gwas_ws$N)
      res_ws$Rb_SumRR <- as.numeric(R_ws %*% res_ws$b_SumRR)
      res_ws$h2_SumRR <- as.numeric(res_ws$b_SumRR %*% R_ws %*% res_ws$b_SumRR)
      ## step2: SumHEM estimates
      Q <- R_ws*median(df_gwas_ws$N) + Matrix::.symDiagonal(ncol(R_ws), x = res_ws$lambda_SumRR)
      res_ws$hv <- res_ws$lambda_SumRR * diag(solve(Q))
      res_ws$hv <- ifelse(res_ws$hv>(1-1e-10), 1-1e-10, res_ws$hv)
      res_ws$sig2b_SumHEM <- res_ws$b_SumRR^2 / (1-res_ws$hv)
      res_ws$sig2b_SumHEM <- ifelse(res_ws$sig2b_SumHEM<1e-10, 1e-10, res_ws$sig2b_SumHEM)
      res_ws$lambda_SumHEM <- res_ws$sig2e / res_ws$sig2b_SumHEM
      res_ws$b_SumHEM <- bigsparser::sp_solve_sym(A = R_SFBM, b = df_gwas_ws$sbeta, add_to_diag = res_ws$lambda_SumHEM/df_gwas_ws$N)
      res_ws$Rb_SumHEM <- as.numeric(R_ws %*% res_ws$b_SumHEM)
      res_ws$h2_SumHEM <- as.numeric(res_ws$b_SumHEM %*% R_ws %*% res_ws$b_SumHEM)
      ## variance of inferred h2 from SumHEM
      H <- res_ws$lambda_SumRR * solve(Q)
      Sig <- res_ws$b_SumRR %*% t(res_ws$b_SumRR) / (1-H)
      Sig[Sig<1e-10 & Sig>0] <- 1e-10
      Sig[Sig>-1e-10 & Sig<0] <- -1e-10
      P <- R_ws %*% Sig
      Mu <- as.numeric(res_ws$b_SumHEM %*% R_ws %*% Sig %*% R_ws %*% res_ws$b_SumHEM)
      res_ws$h2var_SumHEM <- as.numeric(2*sum(diag(P %*% P)) + 4*Mu)
      rm(Q,H,Sig,P,Mu)
      return(as.data.frame(res_ws))
    }
    ## estimate genetic effects
    st <- Sys.time()
    cat(paste0("Estimating genetic effects ...\n"))
    df_est <- NULL
    for (CHR in sort(unique(intersect(df_gwas$chr, df_map$chr)))) {
      st_chr <- Sys.time()
      cat(paste0("Estimating in chromosome ",CHR," ...\t"))
      ## data
      df_gwas_chr <- df_gwas %>% filter(chr==CHR)
      df_map_chr <- df_map %>% filter(chr==CHR) %>% 
        mutate(mapID = rank(pos)) %>% 
        inner_join(df_gwas_chr, by = c("rsid","chr","pos"))
      R_chr <- readRDS(paste0(LD_path,"/LD_chr",CHR,".rds"))
      R_chr@x[is.na(R_chr@x)] <- 1e-10
      ## computation
      df_est_chr <- do.call(rbind, parallel::mclapply(
        1:ceiling(nrow(df_gwas_chr)/ws), estimation_ws,
        ws = ws, df_gwas_chr = df_gwas_chr, df_map_chr = df_map_chr, 
        R_chr = R_chr, h2input = h2input, M = nrow(df_gwas), 
        mc.cores = NCORES
      ))
      df_est_chr <- cbind(df_gwas_chr, df_est_chr)
      df_est <- rbind(df_est, df_est_chr)
      rm(df_gwas_chr,df_map_chr,R_chr,df_est_chr)
      ed_chr <- Sys.time()
      cat(paste0("Time ",difftime(ed_chr,st_chr,units="mins")," min \n"))
    }
    ed <- Sys.time()
    cat(paste0("SumHEM estimation is DONE with Time ",difftime(ed,st,units="mins")," min \n"))
    
    time_end <- date()
    cat(paste0("Analysis finished on ",time_end,"\n"))
    
    return(df_est)
  }

