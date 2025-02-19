
###############################################################################
##### ------------------------- data used to plot ----------------------- #####
###############################################################################

require(dplyr)
workpath = "SumHEM/figures"
##### sim
# simulation results: r2
df_eval_r2_sum <- readRDS(paste0(workpath,"/data/sim/df_eval_r2_sum.rds"))
df_eval_r2_sum_melt <- df_eval_r2_sum %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws")) %>%
  mutate(statistic = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,1],
         method = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,2],
         method = ifelse(method=="true","True",
                         ifelse(method=="gwas","GWAS",
                                ifelse(method=="inf_LDpred","LDpred2-inf",
                                       ifelse(method=="auto_LDpred","LDpred2-auto",
                                              ifelse(method=="1_RRpred","RR",
                                                     ifelse(method=="2_RRpred","SumHEM",NA)))))))
# simulation results: h2
df_eval_h2_sum <- readRDS(paste0(workpath,"/data/sim/df_eval_h2_sum.rds"))
df_eval_h2_sum_melt <- df_eval_h2_sum %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws")) %>%
  mutate(statistic = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,1],
         method = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,2],
         method = ifelse(method=="true","True",
                         ifelse(method=="gwas","GWAS",
                                ifelse(method=="inf_LDpred","LDpred2-inf",
                                       ifelse(method=="auto_LDpred","LDpred2-auto",
                                              ifelse(method=="1_RRpred","RR",
                                                     ifelse(method=="2_RRpred","SumHEM",NA)))))))
# simulation results: effect size cor
df_eval_corabsb_sum <- readRDS(paste0(workpath,"/data/sim/df_eval_corabsb_sum.rds"))
df_eval_corabsb_sum_melt <- df_eval_corabsb_sum %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws")) %>%
  mutate(statistic = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,1],
         method = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,2],
         method = ifelse(method=="true","True",
                         ifelse(method=="gwas","GWAS",
                                ifelse(method=="inf_LDpred","LDpred2-inf",
                                       ifelse(method=="auto_LDpred","LDpred2-auto",
                                              ifelse(method=="1_RRpred","RR",
                                                     ifelse(method=="2_RRpred","SumHEM",NA)))))))
# simulation results: effect direction f1-score
df_eval_f1_sum <- readRDS(paste0(workpath,"/data/sim/df_eval_f1_sum.rds"))
df_eval_f1_sum_melt <- df_eval_f1_sum %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws")) %>%
  mutate(statistic = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,1],
         method = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,2],
         method = ifelse(method=="true","True",
                         ifelse(method=="gwas","GWAS",
                                ifelse(method=="inf_LDpred","LDpred2-inf",
                                       ifelse(method=="auto_LDpred","LDpred2-auto",
                                              ifelse(method=="1_RRpred","RR",
                                                     ifelse(method=="2_RRpred","SumHEM",NA)))))))
# simulation results: beta distribution (mean in bins)
df_eval_meanbins_sum <- readRDS(paste0(workpath,"/data/sim/df_eval_meanbins_sum.rds"))
df_eval_meanbins_sum_melt <- df_eval_meanbins_sum %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws","statistic","method")) %>% 
  mutate(bins = as.integer(stringr::str_remove(variable,"bin_")))
# simulation results: r2 | h2input
df_eval_r2_sum_h2input <- readRDS(paste0(workpath,"/data/sim/df_eval_r2_sum_h2input.rds"))
df_eval_r2_sum_h2input_melt <- df_eval_r2_sum_h2input %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws")) %>%
  mutate(statistic = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,1],
         method = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,2],
         method = ifelse(method=="true","True",
                         ifelse(method=="gwas","GWAS",
                                ifelse(method=="inf_LDpred","LDpred2-inf",
                                       ifelse(method=="auto_LDpred","LDpred2-auto",
                                              ifelse(method=="1_RRpred","RR",
                                                     ifelse(method=="2_RRpred","SumHEM",NA)))))))
# simulation results: h2 | h2input
df_eval_h2_sum_h2input <- readRDS(paste0(workpath,"/data/sim/df_eval_h2_sum_h2input.rds"))
df_eval_h2_sum_h2input_melt <- df_eval_h2_sum_h2input %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws")) %>%
  mutate(statistic = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,1],
         method = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,2],
         method = ifelse(method=="true","True",
                         ifelse(method=="gwas","GWAS",
                                ifelse(method=="inf_LDpred","LDpred2-inf",
                                       ifelse(method=="auto_LDpred","LDpred2-auto",
                                              ifelse(method=="1_RRpred","RR",
                                                     ifelse(method=="2_RRpred","SumHEM",NA)))))))
# simulation results: effect size cor | h2input
df_eval_corabsb_sum_h2input <- readRDS(paste0(workpath,"/data/sim/df_eval_corabsb_sum_h2input.rds"))
df_eval_corabsb_sum_h2input_melt <- df_eval_corabsb_sum_h2input %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws")) %>%
  mutate(statistic = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,1],
         method = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,2],
         method = ifelse(method=="true","True",
                         ifelse(method=="gwas","GWAS",
                                ifelse(method=="inf_LDpred","LDpred2-inf",
                                       ifelse(method=="auto_LDpred","LDpred2-auto",
                                              ifelse(method=="1_RRpred","RR",
                                                     ifelse(method=="2_RRpred","SumHEM",NA)))))))
# simulation results: effect direction f1-score | h2input
df_eval_f1_sum_h2input <- readRDS(paste0(workpath,"/data/sim/df_eval_f1_sum_h2input.rds"))
df_eval_f1_sum_h2input_melt <- df_eval_f1_sum_h2input %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws")) %>%
  mutate(statistic = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,1],
         method = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,2],
         method = ifelse(method=="true","True",
                         ifelse(method=="gwas","GWAS",
                                ifelse(method=="inf_LDpred","LDpred2-inf",
                                       ifelse(method=="auto_LDpred","LDpred2-auto",
                                              ifelse(method=="1_RRpred","RR",
                                                     ifelse(method=="2_RRpred","SumHEM",NA)))))))
# simulation results: beta distribution (mean in bins) | h2input
df_eval_meanbins_sum_h2input <- readRDS(paste0(workpath,"/data/sim/df_eval_meanbins_sum_h2input.rds"))
df_eval_meanbins_sum_h2input_melt <- df_eval_meanbins_sum_h2input %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws","statistic","method")) %>% 
  mutate(bins = as.integer(stringr::str_remove(variable,"bin_")))
# simulation results: r2 | ws
df_eval_r2_sum_ws <- readRDS(paste0(workpath,"/data/sim/df_eval_r2_sum_ws.rds"))
df_eval_r2_sum_ws_melt <- df_eval_r2_sum_ws %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws")) %>%
  mutate(statistic = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,1],
         method = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,2],
         method = ifelse(method=="true","True",
                         ifelse(method=="gwas","GWAS",
                                ifelse(method=="inf_LDpred","LDpred2-inf",
                                       ifelse(method=="auto_LDpred","LDpred2-auto",
                                              ifelse(method=="1_RRpred","RR",
                                                     ifelse(method=="2_RRpred","SumHEM",NA)))))))
# simulation results: h2 | ws
df_eval_h2_sum_ws <- readRDS(paste0(workpath,"/data/sim/df_eval_h2_sum_ws.rds"))
df_eval_h2_sum_ws_melt <- df_eval_h2_sum_ws %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws")) %>%
  mutate(statistic = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,1],
         method = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,2],
         method = ifelse(method=="true","True",
                         ifelse(method=="gwas","GWAS",
                                ifelse(method=="inf_LDpred","LDpred2-inf",
                                       ifelse(method=="auto_LDpred","LDpred2-auto",
                                              ifelse(method=="1_RRpred","RR",
                                                     ifelse(method=="2_RRpred","SumHEM",NA)))))))
# simulation results: effect size cor | ws
df_eval_corabsb_sum_ws <- readRDS(paste0(workpath,"/data/sim/df_eval_corabsb_sum_ws.rds"))
df_eval_corabsb_sum_ws_melt <- df_eval_corabsb_sum_ws %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws")) %>%
  mutate(statistic = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,1],
         method = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,2],
         method = ifelse(method=="true","True",
                         ifelse(method=="gwas","GWAS",
                                ifelse(method=="inf_LDpred","LDpred2-inf",
                                       ifelse(method=="auto_LDpred","LDpred2-auto",
                                              ifelse(method=="1_RRpred","RR",
                                                     ifelse(method=="2_RRpred","SumHEM",NA)))))))
# simulation results: effect direction f1-score | ws
df_eval_f1_sum_ws <- readRDS(paste0(workpath,"/data/sim/df_eval_f1_sum_ws.rds"))
df_eval_f1_sum_ws_melt <- df_eval_f1_sum_ws %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws")) %>%
  mutate(statistic = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,1],
         method = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,2],
         method = ifelse(method=="true","True",
                         ifelse(method=="gwas","GWAS",
                                ifelse(method=="inf_LDpred","LDpred2-inf",
                                       ifelse(method=="auto_LDpred","LDpred2-auto",
                                              ifelse(method=="1_RRpred","RR",
                                                     ifelse(method=="2_RRpred","SumHEM",NA)))))))
# simulation results: beta distribution (mean in bins) | ws
df_eval_meanbins_sum_ws <- readRDS(paste0(workpath,"/data/sim/df_eval_meanbins_sum_ws.rds"))
df_eval_meanbins_sum_ws_melt <- df_eval_meanbins_sum_ws %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws","statistic","method")) %>% 
  mutate(bins = as.integer(stringr::str_remove(variable,"bin_")))
# simulation results: h2 | LDpred2-auto without h2 updating
df_eval_h2_sum_myLDpred <- readRDS(paste0(workpath,"/data/sim/df_eval_h2_sum_myLDpred.rds"))
df_eval_h2_sum_myLDpred_melt <- df_eval_h2_sum_myLDpred %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws")) %>%
  mutate(statistic = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,1],
         method = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,2],
         method = ifelse(method=="true","True",
                         ifelse(method=="gwas","GWAS",
                                ifelse(method=="inf_LDpred","LDpred2-inf",
                                       ifelse(method=="my_LDpred","LDpred2-auto without h2 updating",
                                              ifelse(method=="auto_LDpred","LDpred2-auto with h2 updating",
                                                     ifelse(method=="1_RRpred","RR",
                                                            ifelse(method=="2_RRpred","SumHEM",NA))))))))
# simulation results: estimates of r2
df_eval_eqr2_sum <- readRDS(paste0(workpath,"/data/sim/df_eval_eqr2_sum.rds"))
df_eval_eqr2_sum_melt <- df_eval_eqr2_sum %>% 
  reshape2::melt(id.vars = c("simID","Nsim","N_train","N_test","M","h2","pi","Mcausal","h2input","ws")) %>%
  mutate(statistic = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,1],
         method = stringr::str_split_fixed(variable, pattern = "_", n = 2)[,2],
         method = ifelse(method=="true","True",
                         ifelse(method=="gwas","GWAS",
                                ifelse(method=="inf_LDpred","LDpred2-inf",
                                       ifelse(method=="auto_LDpred","LDpred2-auto",
                                              ifelse(method=="1_RRpred","RR",
                                                     ifelse(method=="2_RRpred","SumHEM",NA)))))))
# simulation results: se of h2_2_RRpred
df_eval_h2se_sum <- readRDS(paste0(workpath,"/data/sim/df_eval_h2se_sum.rds"))
# simulation results: beta
h2h2 = 0.5
pipi = 0.7
simID = 10
df_h2_simID <- readRDS(paste0(workpath,"/data/sim/df_h2_h2",h2h2,"_pi",pipi,"_simID",simID,".rds"))
df_est_simID <- readRDS(paste0(workpath,"/data/sim/df_est_h2",h2h2,"_pi",pipi,"_simID",simID,".rds"))
df_est_simID_melt <- df_est_simID %>% 
  mutate(h2_ldsc = df_h2_simID$h2_ldsc) %>% 
  select(SNP,chr,pos,N,h2_ldsc,b_true,sbeta,b_inf_LDpred,b_auto_LDpred,b_1_RRpred,b_2_RRpred) %>%
  reshape2::melt(id.vars = c("SNP","chr","pos","N","h2_ldsc","b_true")) %>%
  mutate(method = ifelse(variable=="sbeta", "GWAS", 
                         ifelse(variable=="b_inf_LDpred", "LDpred2-inf", 
                                ifelse(variable=="b_auto_LDpred", "LDpred2-auto", 
                                       ifelse(variable=="b_1_RRpred", "RR", 
                                              ifelse(variable=="b_2_RRpred", "SumHEM", NA))))))
##### real
# real analysis results: r2, h2
df_eval_real <- readRDS(paste0(workpath,"/data/real/df_eval.rds"))
df_eval_real_melt <- df_eval_real %>%
  rename(h2_est_inf_LDpred = h2_inf_LDpred, h2_est_auto_LDpred = h2_auto_LDpred, 
         h2_est_1_RRpred = h2_1_RRpred, h2_est_2_RRpred = h2_2_RRpred) %>%
  reshape2::melt(id.vars = c("phenotype","description","source","variable_type",
                             "n_missing","n_non_missing","n_cases","n_controls",
                             "h2_hdl","h2_hdl_se","h2_ldsc","h2_ldsc_se","pi_est")) %>%
  mutate(criteria = stringr::str_split_fixed(variable, pattern = "_", n = 3)[,1], 
         dataset = stringr::str_split_fixed(variable, pattern = "_", n = 3)[,2], 
         dataset = ifelse(dataset=="est", "estimation set", 
                          ifelse(dataset=="val", "validation set", NA)), 
         method = stringr::str_split_fixed(variable, pattern = "_", n = 3)[,3], 
         method = ifelse(method=="inf_LDpred", "LDpred2-inf", 
                         ifelse(method=="auto_LDpred", "LDpred2-auto", 
                                ifelse(method=="1_RRpred", "RR", 
                                       ifelse(method=="2_RRpred", "SumHEM", NA)))))
# real analysis results: beta
phenoID = "20001_1059"
df_est_realID <- readRDS(paste0(workpath,"/data/real/df_est.",phenoID,".rds"))
df_beta_realID <- df_est_realID %>%
  select(SNP,chr,pos) %>%
  group_by(chr) %>%
  summarise(chr_len = length(pos)) %>%
  mutate(chr_tot = cumsum(chr_len)-chr_len) %>%
  right_join(df_est_realID %>% select(SNP,chr,pos,b_gwas,b_gwas_se,sbeta,
                                      b_inf_LDpred,b_auto_LDpred,b_1_RRpred,b_2_RRpred), 
             by = c("chr")) %>%
  group_by(chr) %>%
  mutate(bp = 1:length(pos)) %>%
  ungroup() %>%
  arrange(chr, bp) %>%
  mutate(bp_cum = bp + chr_tot)
df_beta_realID_melt <- df_beta_realID %>%
  mutate(p_gwas = -log10(pchisq((b_gwas/b_gwas_se)^2, df = 1, lower = F))) %>%
  select(SNP,chr,pos,bp_cum,p_gwas,b_inf_LDpred,b_auto_LDpred,b_1_RRpred,b_2_RRpred) %>%
  reshape2::melt(id.vars = c("SNP","chr","pos","bp_cum")) %>%
  mutate(method = ifelse(variable=="p_gwas", "GWAS", 
                         ifelse(variable=="b_inf_LDpred", "LDpred2-inf", 
                                ifelse(variable=="b_auto_LDpred", "LDpred2-auto", 
                                       ifelse(variable=="b_1_RRpred", "RR", 
                                              ifelse(variable=="b_2_RRpred", "SumHEM", NA))))))
df_axis_realID <- df_beta_realID %>%
  group_by(chr) %>%
  summarize(center = (max(bp_cum)+min(bp_cum))/2)

phenoDSCP = df_eval_real$description[df_eval_real$phenotype==phenoID]
pipi = df_eval_real$pi_est[df_eval_real$phenotype==phenoID]
h2h2 = df_eval_real$h2_ldsc[df_eval_real$phenotype==phenoID]

###############################################################################
##### ---------------------------- final plots -------------------------- #####
###############################################################################

##### plot settings
require(dplyr)
require(ggplot2)
workpath = "SumHEM/figures"
ls_color <- list(GWAS = c("#6E348C","#D2BCDE"), 
                 LDpred_auto = c("#276C9E","#A3C9D5"), 
                 LDpred_inf = c("#329845","#AED185"), 
                 LDSC = c("#CCCC99","#4D5D53"), 
                 RRpred_1 = c("#DE7833","#F2BB6B"), 
                 RRpred_2 = c("#C85D4D","#F0B79A"), 
                 True = c("#996600","#CC9966"))
dashed_color = "#4D5D53"
panel_color = "#FAE7D9"

##### ----- ----- ----- ----- figure 1 ----- ----- ----- ----- #####
{
  p_sim_estID <- ggplot() +
    geom_point(data = df_est_simID_melt %>% filter(method!="GWAS") %>% 
                 mutate(method = ifelse(method=="RR","SumRR",method)), 
               mapping = aes(x = b_true, y = value, color = method), 
               alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1, colour = dashed_color, linetype = "dashed") +
    scale_color_manual(values = c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    ggh4x::facet_wrap2(facets = vars(method), nrow = 1, 
                       strip = ggh4x::strip_themed(text_x = ggh4x::elem_list_text(color=c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1])))) +
    scale_x_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = "estimates", x = "true effects", tag = "a") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "white", color = "white"), 
          legend.position = "none", 
          plot.tag = element_text(face = "bold"))
  p_sim_r2_h20.5 <- ggplot() +
    geom_line(data = df_eval_r2_sum_melt %>% filter(statistic=="r2val" & !method%in%c("GWAS","True")) %>% 
                mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                filter(h2==0.5 & pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                group_by(h2,pi,method) %>% summarise(med_value = median(value)), 
              mapping = aes(x = as.factor(pi), y = med_value, group = method, color = method)) +
    geom_errorbar(data = df_eval_r2_sum_melt %>% filter(statistic=="r2val" & !method%in%c("GWAS","True")) %>% 
                    mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                    filter(h2==0.5 & pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                    group_by(h2,pi,method) %>% summarise(med_value = median(value), sd_value = sd(value)), 
                  mapping = aes(x = as.factor(pi), y = med_value, ymin = med_value-sd_value, ymax = med_value+sd_value, color = method), 
                  position = position_dodge(0), width = 0.5) +
    geom_point(data = df_eval_r2_sum_melt %>% filter(statistic=="r2val" & !method%in%c("GWAS","True")) %>% 
                 mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                 filter(h2==0.5 & pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                 group_by(h2,pi,method) %>% summarise(med_value = median(value)), 
               mapping = aes(x = as.factor(pi), y = med_value, color = method)) +
    scale_color_manual(values = c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    scale_x_discrete(labels = c(".1",".3",".5",".7")) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = expression(paste("out-of-sample ",r^2)), x = expression(paste("polygenicity (",pi,")")), color = "methods", tag = "b") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          legend.position = "right", 
          plot.tag = element_text(face = "bold"))
  p_sim_h2_h20.5 <- ggplot() +
    geom_hline(data = data.frame(h2=c(0.5)), mapping = aes(yintercept = h2), colour = dashed_color, linetype = "dashed") +
    geom_line(data = df_eval_h2_sum_melt %>% filter(statistic=="h2" & !method%in%c("GWAS","True")) %>% 
                mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                filter(h2==0.5 & pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                group_by(h2,pi,method) %>% summarise(med_value = median(value)), 
              mapping = aes(x = as.factor(pi), y = med_value, group = method, color = method)) +
    geom_errorbar(data = df_eval_h2_sum_melt %>% filter(statistic=="h2" & !method%in%c("GWAS","True")) %>% 
                    mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                    filter(h2==0.5 & pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                    group_by(h2,pi,method) %>% summarise(med_value = median(value), sd_value = sd(value)), 
                  mapping = aes(x = as.factor(pi), y = med_value, ymin = med_value-sd_value, ymax = med_value+sd_value, color = method), 
                  position = position_dodge(0), width = 0.5) +
    geom_point(data = df_eval_h2_sum_melt %>% filter(statistic=="h2" & !method%in%c("GWAS","True")) %>% 
                 mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                 filter(h2==0.5 & pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                 group_by(h2,pi,method) %>% summarise(med_value = median(value)), 
               mapping = aes(x = as.factor(pi), y = med_value, color = method)) +
    scale_color_manual(values = c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1],ls_color$True[1])) +
    scale_x_discrete(labels = c(".1",".3",".5",".7")) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = expression(paste("inferred ",h^2)), x = expression(paste("polygenicity (",pi,")")), color = "methods", tag = "c") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          legend.position = "right", 
          plot.tag = element_text(face = "bold"))
  p_sim_macrof1_h20.5 <- ggplot() +
    geom_line(data = df_eval_f1_sum_melt %>% filter(statistic=="macrof1" & !method%in%c("GWAS","True")) %>% 
                mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                filter(h2==0.5 & pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                group_by(h2,pi,method) %>% summarise(med_value = median(value)), 
              mapping = aes(x = as.factor(pi), y = med_value, color = method, group = method)) +
    geom_errorbar(data = df_eval_f1_sum_melt %>% filter(statistic=="macrof1" & !method%in%c("GWAS","True")) %>% 
                    mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                    filter(h2==0.5 & pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                    group_by(h2,pi,method) %>% summarise(med_value = median(value), sd_value=sd(value)), 
                  mapping = aes(x = as.factor(pi), y = med_value, color = method, ymin = med_value-sd_value, ymax = med_value+sd_value), 
                  position = position_dodge(0), width = 0.5) +
    geom_point(data = df_eval_f1_sum_melt %>% filter(statistic=="macrof1" & !method%in%c("GWAS","True")) %>% 
                 mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                 filter(h2==0.5 & pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                 group_by(h2,pi,method) %>% summarise(med_value = median(value)), 
               mapping = aes(x = as.factor(pi), y = med_value, color = method)) +
    scale_color_manual(values = c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    scale_x_discrete(labels = c(".1",".3",".5",".7")) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = "macro F1-score", x = expression(paste("polygenicity (",pi,")")), color = "methods", tag = "d") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          legend.position = "right", 
          plot.tag = element_text(face = "bold"))
  p_sim_meanbins_h20.5 <- ggplot(data = df_eval_meanbins_sum_melt %>% 
                                 filter(h2==0.5 & pi==0.7 & !method%in%c("nSNP in bins","True","GWAS")) %>% 
                                   mutate(method = ifelse(method=="RR","SumRR",method))) +
    geom_boxplot(mapping = aes(x = as.factor(bins), y = value, fill = method, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = dashed_color) +
    ggh4x::facet_wrap2(facets = vars(method), nrow = 1,
                       strip = ggh4x::strip_themed(text_x = ggh4x::elem_list_text(color=c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1])))) +
    scale_fill_manual(values = c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    scale_color_manual(values = c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = "mean effects", x = "32 bins", fill =  "methods", tag = "e") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          strip.background = element_rect(fill = "white", color = "white"), 
          legend.position = "none", 
          plot.tag = element_text(face = "bold"))
  p_sim_meanbins_true_h20.5 <- ggplot(data = df_eval_meanbins_sum_melt %>% 
                                  filter(h2==0.5 & pi==0.7 & method%in%c("True"))) +
    geom_boxplot(mapping = aes(x = as.factor(bins), y = value, fill = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = dashed_color) +
    facet_wrap(facets = vars(method), nrow = 1) +
    scale_fill_manual(values = c(ls_color$True[1])) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = "mean effects", x = "32 bins", fill =  "methods", tag = " ") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          axis.text.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.ticks.x = element_blank(), 
          strip.text = element_text(color = "black"), 
          strip.background = element_rect(fill = "white", color = "white"), 
          panel.background = element_rect(fill = panel_color), 
          legend.position = "none", 
          plot.tag = element_text(face = "bold"))
}
p_row1 <- p_sim_estID
p_row2 <- ggpubr::ggarrange(p_sim_r2_h20.5,p_sim_h2_h20.5,p_sim_macrof1_h20.5, ncol = 3, nrow = 1, 
                            common.legend = TRUE, legend = "bottom")
p_row3 <- ggpubr::ggarrange(p_sim_meanbins_h20.5,p_sim_meanbins_true_h20.5, ncol = 2, nrow = 1, widths = c(3.3,1))
figure_1 <- ggpubr::ggarrange(p_row1,p_row2,p_row3, ncol = 1, nrow = 3, heights = c(1.3,1.3,1))
##### ----- ----- ----- ----- figure 2 ----- ----- ----- ----- #####
{
  p_sim_r2_all <- ggplot() +
    geom_line(data = df_eval_r2_sum_melt %>% filter(statistic=="r2val" & !method%in%c("GWAS","True")) %>% 
                mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                filter(pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                group_by(h2,pi,method) %>% summarise(med_value = median(value)) %>% 
                mutate(mylabel = factor(h2, levels = c(0.1,0.5,0.7), 
                                        labels = c(expression(h^2==0.1),expression(h^2==0.5),expression(h^2==0.7)))), 
              mapping = aes(x = as.factor(pi), y = med_value, group = method, color = method)) +
    geom_errorbar(data = df_eval_r2_sum_melt %>% filter(statistic=="r2val" & !method%in%c("GWAS","True")) %>% 
                    mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                    filter(pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                    group_by(h2,pi,method) %>% summarise(med_value = median(value), sd_value = sd(value)) %>% 
                    mutate(mylabel = factor(h2, levels = c(0.1,0.5,0.7), 
                                            labels = c(expression(h^2==0.1),expression(h^2==0.5),expression(h^2==0.7)))), 
                  mapping = aes(x = as.factor(pi), y = med_value, ymin = med_value-sd_value, ymax = med_value+sd_value, color = method), 
                  position = position_dodge(0), width = 0.5) +
    geom_point(data = df_eval_r2_sum_melt %>% filter(statistic=="r2val" & !method%in%c("GWAS","True")) %>% 
                 mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                 filter(pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                 group_by(h2,pi,method) %>% summarise(med_value = median(value)) %>% 
                 mutate(mylabel = factor(h2, levels = c(0.1,0.5,0.7), 
                                         labels = c(expression(h^2==0.1),expression(h^2==0.5),expression(h^2==0.7)))), 
               mapping = aes(x = as.factor(pi), y = med_value, color = method)) +
    scale_color_manual(values = c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    facet_wrap(facets = vars(mylabel), scales = "free_y", labeller = labeller(mylabel=label_parsed)) +
    scale_x_discrete(labels = c(".1",".3",".5",".7")) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = expression(paste("out-of-sample ",r^2)), x = expression(paste("polygenicity (",pi,")")), color = "methods", tag = "a") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          strip.text.x = element_text(size = 10), 
          strip.background = element_rect(fill = "white", color = "white"), 
          legend.position = "right", 
          plot.tag = element_text(face = "bold"))
  p_sim_h2_all <- ggplot() +
    geom_hline(data = data.frame(h2=c(0.1,0.5,0.7)) %>% 
                 mutate(mylabel = factor(h2, levels = c(0.1,0.5,0.7), 
                                         labels = c(expression(h^2==0.1),expression(h^2==0.5),expression(h^2==0.7)))), 
               mapping = aes(yintercept = h2), colour = dashed_color, linetype = "dashed") +
    geom_line(data = df_eval_h2_sum_melt %>% filter(statistic=="h2" & !method%in%c("GWAS","True")) %>% 
                mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                filter(pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                group_by(h2,pi,method) %>% summarise(med_value = median(value)) %>% 
                mutate(mylabel = factor(h2, levels = c(0.1,0.5,0.7), 
                                        labels = c(expression(h^2==0.1),expression(h^2==0.5),expression(h^2==0.7)))), 
              mapping = aes(x = as.factor(pi), y = med_value, group = method, color = method)) +
    geom_errorbar(data = df_eval_h2_sum_melt %>% filter(statistic=="h2" & !method%in%c("GWAS","True")) %>% 
                    mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                    filter(pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                    group_by(h2,pi,method) %>% summarise(med_value = median(value), sd_value = sd(value)) %>% 
                    mutate(mylabel = factor(h2, levels = c(0.1,0.5,0.7), 
                                            labels = c(expression(h^2==0.1),expression(h^2==0.5),expression(h^2==0.7)))), 
                  mapping = aes(x = as.factor(pi), y = med_value, ymin = med_value-sd_value, ymax = med_value+sd_value, color = method), 
                  position = position_dodge(0), width = 0.5) +
    geom_point(data = df_eval_h2_sum_melt %>% filter(statistic=="h2" & !method%in%c("GWAS","True")) %>% 
                 mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                 filter(pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                 group_by(h2,pi,method) %>% summarise(med_value = median(value)) %>% 
                 mutate(mylabel = factor(h2, levels = c(0.1,0.5,0.7), 
                                         labels = c(expression(h^2==0.1),expression(h^2==0.5),expression(h^2==0.7)))), 
               mapping = aes(x = as.factor(pi), y = med_value, color = method)) +
    facet_wrap(facets = vars(mylabel), scales = "free_y", labeller = labeller(mylabel=label_parsed)) +
    scale_color_manual(values = c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1],ls_color$True[1])) +
    scale_x_discrete(labels = c(".1",".3",".5",".7")) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = expression(paste("inferred ",h^2)), x = expression(paste("polygenicity (",pi,")")), color = "methods", tag = "b") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          strip.text.x = element_text(size = 10), 
          strip.background = element_rect(fill = "white", color = "white"), 
          legend.position = "right", 
          plot.tag = element_text(face = "bold"))
  p_sim_macrof1_all <- ggplot() +
    geom_line(data = df_eval_f1_sum_melt %>% filter(statistic=="macrof1" & !method%in%c("GWAS","True")) %>% 
                mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                filter(pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                group_by(h2,pi,method) %>% summarise(med_value = median(value)) %>% 
                mutate(mylabel = factor(h2, levels = c(0.1,0.5,0.7), 
                                        labels = c(expression(h^2==0.1),expression(h^2==0.5),expression(h^2==0.7)))), 
              mapping = aes(x = as.factor(pi), y = med_value, color = method, group = method)) +
    geom_errorbar(data = df_eval_f1_sum_melt %>% filter(statistic=="macrof1" & !method%in%c("GWAS","True")) %>% 
                    mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                    filter(pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                    group_by(h2,pi,method) %>% summarise(med_value = median(value), sd_value=sd(value)) %>% 
                    mutate(mylabel = factor(h2, levels = c(0.1,0.5,0.7), 
                                            labels = c(expression(h^2==0.1),expression(h^2==0.5),expression(h^2==0.7)))), 
                  mapping = aes(x = as.factor(pi), y = med_value, color = method, ymin = med_value-sd_value, ymax = med_value+sd_value), 
                  position = position_dodge(0), width = 0.5) +
    geom_point(data = df_eval_f1_sum_melt %>% filter(statistic=="macrof1" & !method%in%c("GWAS","True")) %>% 
                 mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                 filter(pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                 group_by(h2,pi,method) %>% summarise(med_value = median(value)) %>% 
                 mutate(mylabel = factor(h2, levels = c(0.1,0.5,0.7), 
                                         labels = c(expression(h^2==0.1),expression(h^2==0.5),expression(h^2==0.7)))), 
               mapping = aes(x = as.factor(pi), y = med_value, color = method)) +
    facet_wrap(facets = vars(mylabel), scales = "free_y", labeller = labeller(mylabel=label_parsed)) +
    scale_color_manual(values = c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    scale_x_discrete(labels = c(".1",".3",".5",".7")) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = "macro F1-score", x = expression(paste("polygenicity (",pi,")")), color = "methods", tag = "c") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          strip.text.x = element_text(size = 10), 
          strip.background = element_rect(fill = "white", color = "white"), 
          legend.position = "right", 
          plot.tag = element_text(face = "bold"))
}
figure_2 <- ggpubr::ggarrange(p_sim_r2_all,p_sim_h2_all,p_sim_macrof1_all, ncol = 1, nrow = 3, 
                                common.legend = TRUE, legend = "top")
##### ----- ----- ----- ----- figure 3 ----- ----- ----- ----- #####
{
  p_sim_meanbins_all <- ggplot(data = df_eval_meanbins_sum_melt %>% 
                                 filter(pi==0.7 & !method%in%c("nSNP in bins","True","GWAS")) %>% 
                                 mutate(method = ifelse(method=="RR","SumRR",method)) %>% 
                                 mutate(mylabel = factor(h2, levels = c(0.1,0.5,0.7), 
                                                         labels = c(expression(h^2==0.1),expression(h^2==0.5),expression(h^2==0.7))))) +
    geom_boxplot(mapping = aes(x = as.factor(bins), y = value, fill = method, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = dashed_color) +
    ggh4x::facet_grid2(mylabel ~ method, scales = "free_y", 
                       labeller = labeller(mylabel=label_parsed), 
                       strip = ggh4x::strip_themed(text_x = ggh4x::elem_list_text(color=c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1])))) +
    scale_fill_manual(values = c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    scale_color_manual(values = c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = "mean effects", x = "32 bins", fill =  "methods", tag = "a") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          strip.text.x = element_text(size = 8), 
          # strip.text.y = element_text(size = 10), 
          strip.text.y = element_blank(), 
          strip.background = element_rect(fill = "white", color = "white"), 
          legend.position = "none", 
          plot.tag = element_text(face = "bold"))
  p_sim_meanbins_true <- ggplot(data = df_eval_meanbins_sum_melt %>% 
                                  filter(pi==0.7 & method%in%c("True")) %>% 
                                  mutate(mylabel = factor(h2, levels = c(0.1,0.5,0.7), 
                                                          labels = c(expression(h^2==0.1),expression(h^2==0.5),expression(h^2==0.7))))) +
    geom_boxplot(mapping = aes(x = as.factor(bins), y = value, fill = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = dashed_color) +
    facet_grid(facets = mylabel ~ method, scales = "free_y", labeller = labeller(mylabel=label_parsed)) +
    scale_fill_manual(values = c(ls_color$True[1])) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = "mean effects", x = "32 bins", fill =  "methods", tag = " ") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          axis.text.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.ticks.x = element_blank(), 
          strip.text.x = element_text(size = 8), 
          strip.text.y = element_text(size = 10), 
          strip.background = element_rect(fill = "white", color = "white"), 
          panel.background = element_rect(fill = panel_color), 
          legend.position = "none", 
          plot.tag = element_text(face = "bold"))
}
figure_3 <- ggpubr::ggarrange(p_sim_meanbins_all,p_sim_meanbins_true, ncol = 2, nrow = 1, widths = c(2.5,1))
##### ----- ----- ----- ----- figure 4 ----- ----- ----- ----- #####
{
  p_real_r2 <- ggplot(data = df_eval_real %>% mutate(winner = ifelse(r2_val_auto_LDpred>r2_val_2_RRpred,"LDpred2-auto","SumHEM"))) +
    geom_abline(slope = 1, intercept = 0, color = dashed_color, linetype = "dashed", size = 1, alpha = 0.5) +
    geom_point(mapping = aes(x = r2_val_auto_LDpred, y = r2_val_2_RRpred, colour = winner, size = pi_est)) +
    geom_text(mapping = aes(x = r2_val_auto_LDpred, y = r2_val_2_RRpred, label = description, colour = winner), 
              size = 6, hjust = 0.5, vjust = -1, check_overlap = TRUE) +
    scale_color_manual(values = c(ls_color$LDpred_auto[1], ls_color$RRpred_2[1])) +
    scale_x_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(x = expression(paste("LDpred2-auto ",r^2)), y = expression(paste("SumHEM ",r^2)), 
         color = "Winner", size = expression(paste("Estimated ",pi)), tag = "a") +
    theme_bw() +
    theme(axis.text = element_text(color = "black", size = 18), 
          axis.title = element_text(size = 20), 
          legend.position = "top", 
          legend.title = element_text(size = 18), 
          legend.text = element_text(size = 18), 
          plot.tag = element_text(face = "bold", size = 20))
  p_real_r2_count <- ggplot(data = df_eval_real %>%
                              mutate(winner = ifelse(r2_val_auto_LDpred>r2_val_2_RRpred,"LDpred2-auto","SumHEM"),
                                     pi_group = ifelse(pi_est<0.1,"low","high")) %>%
                              group_by(winner,pi_group) %>%
                              summarise(n_trait = n())) +
    geom_bar(mapping = aes(x = pi_group, y = n_trait, fill = winner),
             stat = "identity", position = "dodge", width = 0.5) +
    scale_fill_manual(values = c(ls_color$LDpred_auto[1], ls_color$RRpred_2[1])) +
    geom_text(mapping = aes(x = pi_group, y = n_trait, group = winner, label = n_trait),
              stat = "identity", position = position_dodge(0.5),
              vjust = -0.5, color = "#4D5D53", size = 3.5, fontface = "bold") +
    lims(y = c(0,155)) +
    labs(y = "# of traits", x = "polygenicity", tag = "b") +
    theme_bw() +
    theme(axis.text = element_text(color = "black", size = 13),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "top",
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 13),
          plot.tag = element_text(face = "bold", size = 20))
  p_real_r2_diff <- ggplot() +
    geom_freqpoly(data = df_eval_real %>% 
                    mutate(r2_diff = r2_val_2_RRpred-r2_val_auto_LDpred, 
                           piclass = ifelse(pi_est>0.1,"high","low")), 
                  mapping = aes(x = r2_diff, color = piclass), 
                  bins = 100, size = 1) +
    geom_vline(xintercept = 0, color = dashed_color, linetype = "dashed") +
    scale_color_manual(values = c(ls_color$RRpred_2[1], ls_color$LDpred_auto[1])) +
    scale_x_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(x = expression(paste("SumHEM ",r^2," - ","LDpred2-auto ",r^2)), y = "# of traits", 
         color = "polygenicity", tag = "c") +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 18), 
          axis.title.x = element_text(size = 18), 
          axis.title.y = element_text(size = 20), 
          legend.position = "top", 
          legend.title = element_text(size = 20), 
          legend.text = element_text(size = 18), 
          plot.tag = element_text(face = "bold", size = 20))
  
  rmse_auto_LDpred = sqrt(mean((df_eval_real$h2_auto_LDpred-df_eval_real$h2_ldsc)^2))
  rmse_2_RRpred = sqrt(mean((df_eval_real$h2_2_RRpred-df_eval_real$h2_ldsc)^2))
  p_real_h2 <- ggplot() +
    geom_text(data = data.frame(method = c("LDpred2-auto","SumHEM"), rmse = c(rmse_auto_LDpred,rmse_2_RRpred)),
              mapping = aes(label = paste0("rmse = ",round(rmse,4))), 
              x = 0, y = 0.6, size = 7, color = "#333333", hjust = 0, vjust = 0) +
    geom_abline(slope = 1, intercept = 0, color = "#4D5D53", linetype = "dashed") +
    geom_point(data = df_eval_real_melt %>% filter(criteria=="h2" & method%in%c("LDpred2-auto","SumHEM")), 
               mapping = aes(x = h2_ldsc, y = value, color = method)) +
    scale_color_manual(values = c(ls_color$LDpred_auto[1],ls_color$RRpred_2[1])) +
    ggh4x::facet_wrap2(facets = vars(method), nrow = 1, 
                       strip = ggh4x::strip_themed(text_x = ggh4x::elem_list_text(color=c(ls_color$LDpred_auto[1],ls_color$RRpred_2[1])))) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(x = expression(paste("LDSC ",h^2)), y = expression(paste("inferred ",h^2)), tag = "d") +
    theme_test() +
    theme(axis.text = element_text(color = "black", size = 16), 
          axis.title = element_text(size = 18), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.line.x.top = element_line(color = "black"), 
          axis.line.y.right = element_line(color = "black"), 
          strip.background = element_rect(fill = "white", color = "white"), 
          strip.text = element_text(size = 18), 
          panel.background = element_blank(), 
          legend.position = "none", 
          plot.tag = element_text(face = "bold", size = 20))
  
}
p_col <- p_real_r2 +
  annotation_custom(grob = ggplotGrob(p_real_r2_count), xmin = 0.18, xmax = 0.27, ymin = 0, ymax = 0.1)
p_row <- ggpubr::ggarrange(p_real_r2_diff,p_real_h2, ncol = 2, nrow = 1)
figure_4 <- ggpubr::ggarrange(p_col,p_row, ncol = 1, nrow = 2, heights = c(10,3))
##### ----- ----- ----- ----- figure 5 ----- ----- ----- ----- #####
{
  p_real_mht <- ggplot() +
    geom_point(data = df_beta_realID_melt %>% filter(method %in% c("GWAS","LDpred2-auto","RR","SumHEM")) %>% 
                 mutate(method = ifelse(method=="RR","SumRR",method)), 
               mapping = aes(x = bp_cum, y = abs(value), color = as.factor(paste0(method,chr%%2)))) +
    scale_color_manual(values = c(ls_color$GWAS[1],ls_color$GWAS[2], 
                                  ls_color$LDpred_auto[1],ls_color$LDpred_auto[2], 
                                  ls_color$RRpred_2[1],ls_color$RRpred_2[2], 
                                  ls_color$RRpred_1[1],ls_color$RRpred_1[2])) +
    scale_x_continuous(label = df_axis_realID$chr, breaks = df_axis_realID$center) +
    ggh4x::facet_wrap2(facets = vars(method), nrow = 4, scales = "free_y", 
                       strip = ggh4x::strip_themed(text_x = ggh4x::elem_list_text(color=c(ls_color$GWAS[1],ls_color$LDpred_auto[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1])))) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(title = expression(paste("self-reported: malignant melanoma (",hat(pi)==0.103,")")), # paste0(phenoDSCP," (pi = ",round(pipi,4),")")
         x = "chromosome", y = "estimates", tag = "a") +
    theme_bw() +
    theme(axis.text = element_text(color = "black", size = 14), 
          axis.title = element_text(size = 18), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          strip.background = element_rect(fill = "white", color = "white"), 
          strip.text = element_text(size = 18), 
          legend.position = "none", 
          plot.title = element_text(face = "bold", size = 18), 
          plot.tag = element_text(face = "bold", size = 20))
  p_real_sca <- ggplot() +
    geom_point(data = df_est_realID, 
               mapping = aes(x = abs(b_1_RRpred), y = abs(b_2_RRpred)), 
               shape = 1) +
    scale_x_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = "SumHEM", x = "SumRR", tag = "b") +
    theme_classic() +
    theme(axis.title = element_text(size = 18), 
          axis.text = element_text(size = 12, color = "black"), 
          axis.text.x = element_text(angle = 30, hjust = 1), 
          legend.position = "none", 
          plot.tag = element_text(face = "bold", size = 20))
  p_real_hv <- ggplot() +
    geom_histogram(data = df_est_realID, mapping = aes(x = hv_RRpred), 
                   color = "black", fill = "white") +
    scale_x_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(y = "# of SNP", x = "hat values", tag = "c") +
    theme_classic() +
    theme(axis.title.x = element_text(size = 18), 
          axis.title.y = element_text(size = 18), 
          axis.text = element_text(size = 12, color = "black"), 
          axis.text.x = element_text(angle = 30, hjust = 1), 
          legend.position = "none", 
          plot.tag = element_text(face = "bold", size = 20))
  rmse_auto_LDpred = sqrt(mean((df_eval_real$h2_auto_LDpred-df_eval_real$h2_ldsc)^2))
  rmse_inf_LDpred = sqrt(mean((df_eval_real$h2_inf_LDpred-df_eval_real$h2_ldsc)^2))
  rmse_1_RRpred = sqrt(mean((df_eval_real$h2_1_RRpred-df_eval_real$h2_ldsc)^2))
  rmse_2_RRpred = sqrt(mean((df_eval_real$h2_2_RRpred-df_eval_real$h2_ldsc)^2))
  p_real_h2 <- ggplot() +
    geom_text(data = data.frame(method = c("LDpred2-auto","LDpred2-inf","SumRR","SumHEM"), rmse = c(rmse_auto_LDpred,rmse_inf_LDpred,rmse_1_RRpred,rmse_2_RRpred)),
              mapping = aes(label = paste0("rmse = ",round(rmse,4))), 
              x = 0, y = 0.6, size = 7, color = "#333333", hjust = 0, vjust = 0) +
    geom_abline(slope = 1, intercept = 0, color = "#4D5D53", linetype = "dashed") +
    geom_point(data = df_eval_real_melt %>% filter(criteria=="h2" & method%in%c("LDpred2-auto","LDpred2-inf","RR","SumHEM")) %>% 
                 mutate(method = ifelse(method=="RR","SumRR",method)), 
               mapping = aes(x = h2_ldsc, y = value, color = method)) +
    scale_color_manual(values = c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    ggh4x::facet_wrap2(facets = vars(method), nrow = 1, 
                       strip = ggh4x::strip_themed(text_x = ggh4x::elem_list_text(color=c(ls_color$LDpred_auto[1],ls_color$LDpred_inf[1],ls_color$RRpred_2[1],ls_color$RRpred_1[1])))) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(x = expression(paste("LDSC ",h^2)), y = expression(paste("inferred ",h^2)), tag = "d") +
    theme_test() +
    theme(axis.text = element_text(color = "black", size = 16), 
          axis.title = element_text(size = 18), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.line.x.top = element_line(color = "black"), 
          axis.line.y.right = element_line(color = "black"), 
          strip.background = element_rect(fill = "white", color = "white"), 
          strip.text = element_text(size = 18), 
          panel.background = element_blank(), 
          legend.position = "none", 
          plot.tag = element_text(face = "bold", size = 20))
}
p_col <- ggpubr::ggarrange(p_real_sca,p_real_hv, ncol = 1, nrow = 2)
p_row <- ggpubr::ggarrange(p_real_mht,p_col, ncol = 2, nrow = 1, widths = c(3,1))
figure_5 <- ggpubr::ggarrange(p_row,p_real_h2, ncol = 1, nrow = 2, heights = c(6,3.5))
##### ----- ----- ----- ----- figure 6 ----- ----- ----- ----- #####
{
  p_sim_r2_h2input <- ggplot() +
    geom_boxplot(data = df_eval_r2_sum_h2input_melt %>% filter(statistic=="r2val" & method%in%c("RR","SumHEM")) %>% 
                   mutate(method = ifelse(method=="RR","SumRR",method)), 
                 mapping = aes(x = as.factor(h2input), y = value, fill = method)) +
    geom_line(data = df_eval_r2_sum_h2input_melt %>% filter(statistic=="r2val" & method%in%c("RR","SumHEM")) %>% 
                mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                group_by(h2,pi,h2input,method) %>% summarise(med_value = median(value)), 
              mapping = aes(x = as.factor(h2input), y = med_value, group = method, color = method)) +
    geom_point(data = df_eval_r2_sum_h2input_melt %>% filter(statistic=="r2val" & method%in%c("RR","SumHEM")) %>% 
                 mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                 group_by(h2,pi,h2input,method) %>% summarise(med_value = median(value)), 
               mapping = aes(x = as.factor(h2input), y = med_value, color = method)) +
    scale_fill_manual(values = c(ls_color$RRpred_2[2],ls_color$RRpred_1[2])) +
    scale_color_manual(values = c(ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    scale_x_discrete(labels = c(".3",".4",".5",".6",".7")) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = expression(paste("out-of-sample ",r^2)), x = expression(paste("input of ",h^2)), color = NULL, fill = NULL, tag = "a") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          legend.position = "right", 
          plot.tag = element_text(face = "bold"))
  p_sim_h2_h2input <- ggplot() +
    geom_hline(data = data.frame(h2=c(0.5)), mapping = aes(yintercept = h2), colour = dashed_color, linetype = "dashed") +
    geom_boxplot(data = df_eval_h2_sum_h2input_melt %>% filter(statistic=="h2" & method%in%c("RR","SumHEM")) %>% 
                   mutate(method = ifelse(method=="RR","SumRR",method)), 
                 mapping = aes(x = as.factor(h2input), y = value, fill = method)) +
    geom_line(data = df_eval_h2_sum_h2input_melt %>% filter(statistic=="h2" & method%in%c("RR","SumHEM")) %>% 
                mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                group_by(h2,pi,h2input,method) %>% summarise(med_value = median(value)), 
              mapping = aes(x = as.factor(h2input), y = med_value, group = method, color = method)) +
    geom_point(data = df_eval_h2_sum_h2input_melt %>% filter(statistic=="h2" & method%in%c("RR","SumHEM")) %>% 
                 mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                 group_by(h2,pi,h2input,method) %>% summarise(med_value = median(value)), 
               mapping = aes(x = as.factor(h2input), y = med_value, color = method)) +
    scale_fill_manual(values = c(ls_color$RRpred_2[2],ls_color$RRpred_1[2])) +
    scale_color_manual(values = c(ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    scale_x_discrete(labels = c(".3",".4",".5",".6",".7")) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = expression(paste("inferred ",h^2)), x = expression(paste("input of ",h^2)), color = NULL, fill = NULL, tag = "b") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          legend.position = "right", 
          plot.tag = element_text(face = "bold"))
  p_sim_macrof1_h2input <- ggplot() +
    geom_boxplot(data = df_eval_f1_sum_h2input_melt %>% filter(statistic=="macrof1" & method%in%c("RR","SumHEM")) %>% 
                   mutate(method = ifelse(method=="RR","SumRR",method)), 
                 mapping = aes(x = as.factor(h2input), y = value, fill = method)) +
    geom_line(data = df_eval_f1_sum_h2input_melt %>% filter(statistic=="macrof1" & method%in%c("RR","SumHEM")) %>% 
                mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                group_by(h2,pi,h2input,method) %>% summarise(med_value = median(value)), 
              mapping = aes(x = as.factor(h2input), y = med_value, group = method, color = method)) +
    geom_point(data = df_eval_f1_sum_h2input_melt %>% filter(statistic=="macrof1" & method%in%c("RR","SumHEM")) %>% 
                 mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                 group_by(h2,pi,h2input,method) %>% summarise(med_value = median(value)), 
               mapping = aes(x = as.factor(h2input), y = med_value, color = method)) +
    scale_fill_manual(values = c(ls_color$RRpred_2[2],ls_color$RRpred_1[2])) +
    scale_color_manual(values = c(ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    scale_x_discrete(labels = c(".3",".4",".5",".6",".7")) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = "macro F1-score", x = expression(paste("input of ",h^2)), color = NULL, fill = NULL, tag = "c") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          legend.position = "right", 
          plot.tag = element_text(face = "bold"))
  p_sim_meanbins_h2input <- ggplot(data = df_eval_meanbins_sum_h2input_melt %>% 
                                     filter(method%in%c("RR","SumHEM")) %>% 
                                     mutate(method = ifelse(method=="RR","SumRR",method))) +
    geom_boxplot(mapping = aes(x = as.factor(bins), y = value, fill = method, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "#4D5D53") +
    ggh4x::facet_grid2(paste0("input = ",h2input) ~ method, scales = "fixed",
                       strip = ggh4x::strip_themed(text_x = ggh4x::elem_list_text(color=c(ls_color$RRpred_2[1],ls_color$RRpred_1[1])))) +
    scale_fill_manual(values = c(ls_color$RRpred_2[2],ls_color$RRpred_1[2])) +
    scale_color_manual(values = c(ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = "mean effects", x = "32 bins", fill = "methods", tag = "d") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          strip.text.y = element_blank(), 
          strip.background = element_rect(fill = "white", color = "white"), 
          legend.position = "none", 
          plot.tag = element_text(face = "bold"))
  p_sim_meanbins_true_h2input <- ggplot(data = df_eval_meanbins_sum_h2input_melt %>% 
                                          filter(pi==0.7 & method%in%c("True"))) +
    geom_boxplot(mapping = aes(x = as.factor(bins), y = value, fill = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "#4D5D53") +
    facet_grid(facets = paste0("input = ",h2input) ~ method, scales = "fixed") +
    scale_fill_manual(values = c(ls_color$True[1])) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = "mean effects", x = "32 bins", fill = "methods", tag = " ") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          axis.text.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.ticks.x = element_blank(), 
          strip.text = element_text(color = "black"), 
          strip.background = element_rect(fill = "white", color = "white"), 
          panel.background = element_rect(fill = panel_color), 
          legend.position = "none", 
          plot.tag = element_text(face = "bold"))
}
p_col <- ggpubr::ggarrange(p_sim_r2_h2input,p_sim_h2_h2input,p_sim_macrof1_h2input, ncol = 1, nrow = 3, 
                           common.legend = TRUE, legend = "bottom")
p_row <- ggpubr::ggarrange(p_sim_meanbins_h2input,p_sim_meanbins_true_h2input, ncol = 2, nrow = 1, widths = c(1.65,1))
figure_6 <- ggpubr::ggarrange(p_col,p_row, ncol = 2, nrow = 1, widths = c(1.2,3))
##### ----- ----- ----- ----- figure 7 ----- ----- ----- ----- #####
{
  p_sim_r2_ws <- ggplot() +
    geom_boxplot(data = df_eval_r2_sum_ws_melt %>% filter(statistic=="r2val" & method%in%c("RR","SumHEM")) %>% 
                   mutate(method = ifelse(method=="RR","SumRR",method)), 
                 mapping = aes(x = as.factor(ws), y = value, fill = method)) +
    geom_line(data = df_eval_r2_sum_ws_melt %>% filter(statistic=="r2val" & method%in%c("RR","SumHEM")) %>% 
                mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                group_by(h2,pi,h2input,ws,method) %>% summarise(med_value = median(value)), 
              mapping = aes(x = as.factor(ws), y = med_value, group = method, color = method)) +
    geom_point(data = df_eval_r2_sum_ws_melt %>% filter(statistic=="r2val" & method%in%c("RR","SumHEM")) %>% 
                 mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                 group_by(h2,pi,h2input,ws,method) %>% summarise(med_value = median(value)), 
               mapping = aes(x = as.factor(ws), y = med_value, color = method)) +
    scale_fill_manual(values = c(ls_color$RRpred_2[2],ls_color$RRpred_1[2])) +
    scale_color_manual(values = c(ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = expression(paste("out-of-sample ",r^2)), x = "window size", color = NULL, fill = NULL, tag = "a") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          legend.position = "right", 
          plot.tag = element_text(face = "bold"))
  p_sim_h2_ws <- ggplot() +
    geom_hline(data = data.frame(h2=c(0.5)), mapping = aes(yintercept = h2), colour = "#4D5D53", linetype = "dashed") +
    geom_boxplot(data = df_eval_h2_sum_ws_melt %>% filter(statistic=="h2" & method%in%c("RR","SumHEM") & value>0) %>% 
                   mutate(method = ifelse(method=="RR","SumRR",method)), 
                 mapping = aes(x = as.factor(ws), y = value, fill = method)) +
    geom_line(data = df_eval_h2_sum_ws_melt %>% filter(statistic=="h2" & method%in%c("RR","SumHEM") & value>0) %>% 
                mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                group_by(h2,pi,h2input,ws,method) %>% summarise(med_value = median(value)), 
              mapping = aes(x = as.factor(ws), y = med_value, group = method, color = method)) +
    geom_point(data = df_eval_h2_sum_ws_melt %>% filter(statistic=="h2" & method%in%c("RR","SumHEM") & value>0) %>% 
                 mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                 group_by(h2,pi,h2input,ws,method) %>% summarise(med_value = median(value)), 
               mapping = aes(x = as.factor(ws), y = med_value, color = method)) +
    scale_fill_manual(values = c(ls_color$RRpred_2[2],ls_color$RRpred_1[2])) +
    scale_color_manual(values = c(ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = expression(paste("inferred ",h^2)), x = "window size", color = NULL, fill = NULL, tag = "b") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          legend.position = "right", 
          plot.tag = element_text(face = "bold"))
  p_sim_macrof1_ws <- ggplot() +
    geom_boxplot(data = df_eval_f1_sum_ws_melt %>% filter(statistic=="macrof1" & method%in%c("RR","SumHEM")) %>% 
                   mutate(method = ifelse(method=="RR","SumRR",method)), 
                 mapping = aes(x = as.factor(ws), y = value, fill = method)) +
    geom_line(data = df_eval_f1_sum_ws_melt %>% filter(statistic=="macrof1" & method%in%c("RR","SumHEM")) %>% 
                mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                group_by(h2,pi,h2input,ws,method) %>% summarise(med_value = median(value)), 
              mapping = aes(x = as.factor(ws), y = med_value, group = method, color = method)) +
    geom_point(data = df_eval_f1_sum_ws_melt %>% filter(statistic=="macrof1" & method%in%c("RR","SumHEM")) %>% 
                 mutate(method = ifelse(method=="RR","SumRR",method)) %>%
                 group_by(h2,pi,h2input,ws,method) %>% summarise(med_value = median(value)), 
               mapping = aes(x = as.factor(ws), y = med_value, color = method)) +
    scale_fill_manual(values = c(ls_color$RRpred_2[2],ls_color$RRpred_1[2])) +
    scale_color_manual(values = c(ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = "macro F1-score", x = "window size", color = NULL, fill = NULL, tag = "c") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          legend.position = "right", 
          plot.tag = element_text(face = "bold"))
  p_sim_meanbins_ws <- ggplot(data = df_eval_meanbins_sum_ws_melt %>% 
                                filter(method%in%c("RR","SumHEM")) %>% 
                                mutate(method = ifelse(method=="RR","SumRR",method))) +
    geom_boxplot(mapping = aes(x = as.factor(bins), y = value, fill = method, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "#4D5D53") +
    ggh4x::facet_grid2(paste0("window size = ",ws) ~ method, scales = "fixed",
                       strip = ggh4x::strip_themed(text_x = ggh4x::elem_list_text(color=c(ls_color$RRpred_2[1],ls_color$RRpred_1[1])))) +
    scale_fill_manual(values = c(ls_color$RRpred_2[2],ls_color$RRpred_1[2])) +
    scale_color_manual(values = c(ls_color$RRpred_2[1],ls_color$RRpred_1[1])) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = "mean effects", x = "32 bins", fill =  "methods", tag = "d") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          strip.text.y = element_blank(), 
          strip.background = element_rect(fill = "white", color = "white"), 
          legend.position = "none", 
          plot.tag = element_text(face = "bold"))
  p_sim_meanbins_true_ws <- ggplot(data = df_eval_meanbins_sum_ws_melt %>% 
                                     filter(pi==0.7 & method%in%c("True"))) +
    geom_boxplot(mapping = aes(x = as.factor(bins), y = value, fill = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = dashed_color) +
    facet_grid(facets = paste0("window size = ",ws) ~ method, scales = "fixed") +
    scale_fill_manual(values = c(ls_color$True[1])) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = "mean effects", x = "32 bins", fill =  "methods", tag = " ") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          axis.text.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.ticks.x = element_blank(), 
          strip.text = element_text(color = "black"), 
          strip.background = element_rect(fill = "white", color = "white"), 
          panel.background = element_rect(fill = panel_color), 
          legend.position = "none", 
          plot.tag = element_text(face = "bold"))
}
p_col <- ggpubr::ggarrange(p_sim_r2_ws,p_sim_h2_ws,p_sim_macrof1_ws, ncol = 1, nrow = 3, 
                           common.legend = TRUE, legend = "bottom")
p_row <- ggpubr::ggarrange(p_sim_meanbins_ws,p_sim_meanbins_true_ws, ncol = 2, nrow = 1, widths = c(1.65,1))
figure_7 <- ggpubr::ggarrange(p_col,p_row, ncol = 2, nrow = 1, widths = c(1.2,3))
##### ----- ----- ----- ----- figure 8 ----- ----- ----- ----- #####
{
  p_sim_h2se <- ggplot() +
    geom_boxplot(data = df_eval_h2se_sum %>% filter(pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                   mutate(mylabel = factor(h2, levels = c(0.1,0.5,0.7), 
                                           labels = c(expression(h^2==0.1),expression(h^2==0.5),expression(h^2==0.7)))), 
                 mapping = aes(x = as.factor(pi), y = h2se_2_RRpred), fill = ls_color$RRpred_2[2]) +
    geom_line(data = df_eval_h2_sum %>% filter(pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                group_by(h2,pi) %>% summarise(h2sd_2_RRpred = sd(h2_2_RRpred)) %>% 
                mutate(mylabel = factor(h2, levels = c(0.1,0.5,0.7), 
                                        labels = c(expression(h^2==0.1),expression(h^2==0.5),expression(h^2==0.7)))), 
              mapping = aes(x = as.factor(pi), y = h2sd_2_RRpred, group = h2), 
              color = ls_color$RRpred_2[1]) +
    geom_point(data = df_eval_h2_sum %>% filter(pi%in%c(0.1,0.3,0.5,0.7)) %>% 
                 group_by(h2,pi) %>% summarise(h2sd_2_RRpred = sd(h2_2_RRpred)) %>% 
                 mutate(mylabel = factor(h2, levels = c(0.1,0.5,0.7), 
                                         labels = c(expression(h^2==0.1),expression(h^2==0.5),expression(h^2==0.7)))), 
               mapping = aes(x = as.factor(pi), y = h2sd_2_RRpred), 
               color = ls_color$RRpred_2[1]) +
    facet_wrap(facets = vars(mylabel), scales = "free_y", labeller = labeller(mylabel=label_parsed)) +
    scale_x_discrete(labels = c(".1",".3",".5",".7")) +
    scale_y_continuous(labels = function(x) gsub('^0\\.', '\\.', as.character(x))) +
    labs(y = expression(paste("se of inferred ",h^2)), x = expression(paste("polygenicity (",pi,")")), tag = "a") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          strip.background = element_rect(fill = "white", color = "white"), 
          legend.position = "none", 
          plot.tag = element_text(face = "bold"))
  p_sim_updateh2 <- ggplot(data = df_eval_h2_sum_myLDpred_melt %>% 
                             filter(method%in%c("SumHEM","LDpred2-auto without h2 updating","LDpred2-auto with h2 updating")) %>% 
                             mutate(mylabel = factor(pi, levels = c(0.1,0.3,0.5,0.7), 
                                                     labels = c(expression(pi==0.1),expression(pi==0.3),expression(pi==0.5),expression(pi==0.7))))) +
    geom_hline(mapping = aes(yintercept = h2), colour = dashed_color, linetype = "dashed") +
    geom_boxplot(mapping = aes(x = method, y = value, color = method)) +
    facet_wrap(facets = vars(mylabel), nrow = 1, labeller = labeller(mylabel=label_parsed)) +
    labs(y = expression(paste("inferred ",h^2)), x = "method", tag = "b") +
    theme_classic() +
    scale_color_manual(values = c(ls_color$LDpred_auto[1], ls_color$LDpred_auto[2], ls_color$RRpred_2[1])) +
    scale_fill_manual(values = c(ls_color$LDpred_auto[1], ls_color$LDpred_auto[2], ls_color$RRpred_2[1])) +
    theme(axis.text = element_text(color = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text.x = element_text(size = 10),
          strip.background = element_rect(fill = "white", color = "white"), 
          plot.tag = element_text(face = "bold"))
}
figure_8 <- ggpubr::ggarrange(p_sim_h2se,p_sim_updateh2, ncol = 1)












