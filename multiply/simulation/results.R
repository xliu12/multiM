
library(tidyverse)
options(dplyr.print_max = 100, pillar.print_max = 50, dplyr.width = Inf)
library(glue)
library(cowplot)
library(ggpattern)
library(GGally)
# library(latex2exp)


# # removing files
# rmFile <- grep(".Rout", dir("RData"),value = T) %>% paste0("RData/",.)
# file.remove(rmFile)
# 
# # renaming/moving files
# fromFile <- grep("job_[0-9]+.RData", dir("RData"),value = T) %>% paste0("RData/",.)
# toFile <- grep("job_[0-9]+.RData", dir("RData"),value = T) %>% paste0("RData/num_c5/",.)

# file.rename(fromFile, toFile)

# scp file transfer --------
# scp -r xliu19@ls6.tacc.utexas.edu:/work/08878/xliu19/ls6/multiply/RData/ .

# Nonrandomized ---------------------------


# Reps 501 to 1000 ---- 

dir("RData")

dirpath <- "RData/FEglm.RData"
dirpath <- "RData/linear_Fit_FEglm_SLglm.RData"
dirpath <- "RData/linear_Fit_ySLglm.RData"


load(dirpath)

# dirpath_s <- unlist(mget(paste0("dirpath", 1:2)))

dirpath_s <- dir("RData") %>% 
  grep("xz.*_linear_Fit_glm_icc.*_rep1_1000.*", ., value = TRUE)
  # grep("job_B[qrqt]?[qm]?[qy]?[qzero]?[qthree]?", ., value = T) %>% grep("20c",., value = T) %>% paste0("RData/",.) %>% 
  # c(., paste0("RData/reps501to1000/", dir("RData/reps501to1000/")) )



# ...........------------

jobconds <- 1:nrow(condition)


dirpath_s <- dirpath

dflist <-vector(mode = "list", length = nrow(condition))


# (skip) Results reading -------------------

for(i in 1:nrow(condition)) {
  res_i <- vector(mode = "list", length = length(dirpath_s))
  
  for(d in 1:length(dirpath_s)) {
    
    dirpath <- dirpath_s[d]
    if (dir.exists(dirpath)) {
      res_meD <- vector(mode = "list", length = length(dir(dirpath)))
      length_dirpath <- length(dir(dirpath)) 
    }
    if ( (!dir.exists(dirpath)) & (file.exists(dirpath)) ) {
      res_meD <- vector(mode = "list", length = length((dirpath)))
      length_dirpath <- length((dirpath)) 
    }
    
    
    for(f in 1:length_dirpath ) {
      if (dir.exists(dirpath) & file.exists(paste0(dirpath, dir(dirpath)[f])) ) {
        load(paste0(dirpath, dir(dirpath)[f]))
      }
      if ( file.exists(dirpath) ) {
        load(dirpath)
      }
      if(length(reslist) >= i) {
        if( !is.null(reslist[[i]]) ) {
          
          res <- reslist[[i]]
          estimates <- do.call(rbind, res)
          # estimates <- do.call(rbind, lapply(res, `[[`, "results"))
          
          res_meD[[f]] <- estimates #simu_cond_i
        }
      }
    }
    res_di <- do.call(rbind, res_meD)
    
    if(length(res_di) > 1) {
      res_i[[d]] <- data.frame(cond = i, res_di, row.names = NULL)
    }
  }
  
  dflist[[i]] <- res_i
}

dflst <- vector(mode = "list", length = nrow(condition))
for(i in 1:nrow(condition)) {
  dflst[[i]] <- do.call(rbind, dflist[[i]])
}

resdf <- do.call(rbind, dflst)
# .............. ----------------------

library(readr)
dirpath <- dirpath_s[1]

read.one <- function(dirpath) {
  load(file.path("RData", dirpath))
  # sapply(reslist, length) > 0
  res_bind <- lapply(
    reslist[sapply(reslist, length) > 0], 
    function(l){
      do.call(rbind, l[sapply(l, length) > 1])
    })
  
  resdf <- do.call(rbind, res_bind)
  resdf1 <- resdf
  if (!"x_z" %in% colnames(resdf)) {
    resdf$x_z <- 0
    resdf1 <- resdf %>% select(colnames(condition), everything())
  }
  return(resdf1)
}
resdf <- lapply(dirpath_s, read.one)
resdf1 <- do.call(rbind, resdf)
# save(resdf1, file = "./resdf1.RData")

library(tidyverse)
options(dplyr.print_max = 100, pillar.print_max = 50, dplyr.width = Inf)

# meD ----------------

res_df <- resdf1 %>%
  group_by(num_clust, clust_size, quadratic.A, quadratic.M, quadratic.Y, icc, Yfamily, x_z,
           cluster_a, cluster_m, cluster_y, 
           Fit,
           effname) %>% 
  mutate(
    trueval = mean(true_val) 
    ) %>% 
  mutate(
    rbias = across(contains("est_"), ~. / trueval - 1),
    bias = across(contains("est_"), ~. - trueval),
    # rbias_qlogis = across(contains("est_"), ~(qlogis(.) / qlogis(trueval) - 1)),
    # bias_qlogis = across(contains("est_"), ~log(./(1-.)) - log(trueval/(1-trueval))),
    MSE = across(contains("est_"), ~(. - trueval)^2),
    nreps = n(),
    multi_cover = (ci1_multi < (trueval)) & (ci2_multi > (trueval))
    ) 

# save(res_df, file = "simulation/res_df.RData")

res_df <- do.call(data.frame, res_df)


res_summ <- res_df %>%
  group_by(num_clust, clust_size, quadratic.A, quadratic.M, quadratic.Y, icc, Yfamily, x_z,
           cluster_a, cluster_m, cluster_y, 
           Fit,
           effname) %>% 
  summarise( across( c("nreps", contains("trueval"), contains("bias"), contains("MSE"), contains("cover")), 
                    ~mean(., na.rm=T) )
             )  %>%
  select(group_cols(), effname, contains("rbias"), contains("bias"), contains("cover"), contains("MSE"), everything()) %>% 
  mutate(effname = factor(effname))






res_summ <- res_summ %>% mutate(
  across(contains("MSE"), .fns = list(sqrt = ~sqrt(.)
                                      # , sqrtN = ~sqrt(. * n)
                                      ), 
         .names ="{.fn}{.col}" )
)
# _reg aggregate over a_model (because for each combination of methods, there may be data generated by different seeds in parallel)
res_summ1 <- res_summ %>% 
  group_by(num_clust, clust_size, quadratic.A, quadratic.M, quadratic.Y, icc, Yfamily, x_z,
           cluster_m, cluster_y, # cluster_a, 
           Fit,
           effname) %>% 
  mutate(
    across(ends_with("_reg"), ~mean(.))
  )
# _rmpw aggregate over y_model 
res_summ1 <- res_summ1 %>% 
  ungroup() %>% 
  group_by(num_clust, clust_size, quadratic.A, quadratic.M, quadratic.Y, icc, Yfamily, x_z,
           cluster_m,  cluster_a, #cluster_y, #
           Fit,
           effname) %>% 
  mutate(
    across(ends_with("_rmpw"), ~mean(.))
  )

res_long <- res_summ1 %>% select(!ends_with("_cover")) %>%
  pivot_longer(
  # cols = c(ends_with("_rbias"), ends_with("_cover"), ends_with("_MSE"), ends_with("_bias")), 
  cols = c(ends_with("_multi"), ends_with("_reg"), ends_with("_rmpw")), # not the weighting estimator
  names_to = c("performance", "estimator"),
  names_pattern = "(.*).est_(.*)"
  # names_sep = "_"
)




# save .csv -----------------
write_csv(res_summ, file = "Tables_Figs/res_summ.csv")
write_csv(res_summ1, file = "Tables_Figs/res_summ1.csv")

write_csv(res_long, file = "Tables_Figs/res_long.csv")


save(res_summ, file = "simulation/res_summ.RData")



# Load .csv ----
library(readr)
res_summ <- read_csv("Tables_Figs/res_summ.csv")
res_long <- read_csv("Tables_Figs/res_long.csv")

# Re-labelling --------
res_long1 <- res_long %>% mutate(
  # true_model = case_when(
  #   quadratic.R+quadratic.tt+quadratic.M+quadratic.Y==0 ~ "all linear",
  #   quadratic.R+quadratic.tt+quadratic.M+quadratic.Y==4 ~ "all quadratic",
  #   ((quadratic.R+quadratic.tt) %in% c(0) )&
  #     ((quadratic.M) %in% c(1) ) &
  #     ((quadratic.Y) %in% c(0) ) ~ "only M quadratic",
  #   ((quadratic.R+quadratic.tt) %in% c(0) )&
  #     ((quadratic.M) %in% c(0) ) &
  #     ((quadratic.Y) %in% c(1) ) ~ "only Y quadratic",
  #   ((quadratic.R+quadratic.tt) %in% c(2) )&
  #     ((quadratic.M) %in% c(0) ) &
  #     ((quadratic.Y) %in% c(0) ) ~ "only R, T quadratic",
  #   TRUE ~ "other"
  # ),
  # true_model = factor(true_model, levels = c("all linear", "only R, T quadratic", "only M quadratic", "only Y quadratic", "all quadratic", "other")),
  performance = factor(performance),
  sample_size = factor(num_clust * clust_size),
  J = factor(num_clust), nj = factor(clust_size),
  # Estimator = factor(case_when(
  #   str_detect(estimator, "tml") ~ "TML",
  #   str_detect(estimator, "onestep") ~ "One-step"
  # )),
  # Nuisance_estimation = factor(case_when(
  #   Fit=="linear" ~ "parametric",
  #   Fit=="mlr" ~ "machine"
  # ), levels = c("parametric", "machine")),
  cluster_a = case_match(cluster_a,
                         "noncluster" ~ "SL", .default = cluster_a),
  cluster_m = case_match(cluster_m,
                         "noncluster" ~ "SL", .default = cluster_m),
  cluster_y = case_match(cluster_y,
                         "noncluster" ~ "SL", .default = cluster_y),
  A_mod = ifelse(cluster_a %in% c("FE","RE"), "aCluster", "aSL"),
  M_mod = ifelse(cluster_m %in% c("FE","RE"), "mCluster", "mSL"),
  Y_mod = ifelse(cluster_y %in% c("FE","RE"), "yCluster", "ySL"),
  abs.value = abs(value),
  abs.performance = case_when(
    as.character(performance)=="bias" ~ "|Bias|",
    as.character(performance)=="rbias" ~ "|Rel.Bias|",
    as.character(performance)=="sqrtMSE" ~ "RMSE",
    TRUE ~ as.character(performance)
  )
) %>% 
  mutate(
    abs.performance = factor(abs.performance),
    Estimator = case_match(estimator, 
                           "reg" ~ "regression", "rmpw" ~ "weighting", "multi" ~ "multiply-robust", 
                           .ptype = factor(levels = c("multiply-robust", "regression", "weighting")))
    ) %>% 
  mutate(A_fit = case_match(cluster_a, 
                            "FE" ~ "Afixed", "RE" ~ "Arand", "SL" ~ "Asl"), 
         M_fit = case_match(cluster_m, 
                            "FE" ~ "Mfixed", "RE" ~ "Mrand", "SL" ~ "Msl"), 
         Y_fit = case_match(cluster_y, 
                            "FE" ~ "Yfixed", "RE" ~ "Yrand", "SL" ~ "Ysl")
         )

write_csv(res_long1, file = "Tables_Figs/res_long1.csv")

# Table ----------------------------

# across sample sizes 
aggtab <- res_long1 %>% 
  filter(str_detect(effname, "Y\\(") ) %>% 
  group_by(
  J, nj, 
  cluster_a, cluster_m, cluster_y,
  A_fit, M_fit, Y_fit,
  A_mod, M_mod, Y_mod,
  Estimator, estimator,
  # effname, 
  abs.performance
) %>% summarise(
  across(abs.value, list(mean = ~mean(.), min = ~min(.), max = ~max(.)))
) %>% 
  filter(abs.performance != "MSE" )

aggtab_wide <- aggtab %>% pivot_wider(
  names_from = c(abs.performance),
  # names_from = c(abs.performance, estimator),
  values_from = contains("abs.value"),
  # names_prefix = c("n_"),
  # names_glue = "{abs.performance}_J{J}_nj{nj}" # "n{n}_{abs.performance}_{.value}" 
  names_vary = "fastest"
  ,names_sort = T
) #%>% mutate(across(starts_with("MSE_"), ~(.*1000), .names ="1000{.col}") )

aggtab_wide <- aggtab_wide %>% 
  select(# effname, 
         J, nj, Estimator, #A_mod, M_mod, Y_mod, 
         cluster_a, cluster_m, cluster_y,
         contains("Rel.Bias"), contains("RMSE"), everything()) %>% 
  arrange(Estimator)

write_csv(aggtab, file = "Tables_Figs/aggtab.csv")

write_csv(aggtab_wide, file = "Tables_Figs/aggtab_wide.csv")




# save.image here ---------------------
save.image("simulation/results.RData")


# PLOT ----------------
load("simulation/results.RData")

# theme figure -----------
fig.theme <- 
  theme_bw() +
  theme(panel.grid.minor.y = element_line(size = 0),
        panel.grid.major.x = element_line(size = 0),
        panel.grid.major.y = element_line(size = 0.2, lineend = "round", color = "grey", linetype = "solid"),
        # axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.position = "top",
        legend.spacing.x = unit(0.1, "mm"),
        legend.key.height = unit(15, "mm"),
        legend.key.size = unit(10, "mm"))



##  cluster vs SL ----
plotdf <- res_long1 %>% 
  filter(
    abs.performance %in% c("|Rel.Bias|", "RMSE"), 
    str_detect(effname, "Y\\(")
  # effname %in% c("Y(1,gm1(1),gm2(0))") # c( "Y(0,gmjo(0))" ), # (!true_model %in% c("other"))
  , J == 100, nj == 10 #, cluster_m %in% c("FE", "RE")
  )
# plotdf <- aggtab %>% 
  # filter(str_detect(abs.performance, "R"), 
  #        nj == 10, J == 100
  #        # str_detect(effname, "Y\\(")
  #        )

fig.rbias <- ggplot(plotdf) +
  ggtitle( glue("{unique(plotdf$J)} clusters, cluster size {unique(plotdf$nj)}") ) + #, ICC = {unique(plotdf$icc)}, association of X and U = {unique(plotdf$x_z)};  estimand: {unique(plotdf$effname)}
  facet_grid(Estimator ~abs.performance , 
                #A_fit ~ Y_fit, #Y_true + Clus_true + Trt_true
             labeller = labeller(
               # cluster_a = label_both, cluster_y = label_both , 
               abs.performance = label_value, Estimator = label_both
                                 ), 
             scales = "free") +
  # geom_point(
  geom_boxplot(
  aes( x = abs.value, 
    # xmin = , xmax = max(abs.value), xlower=quantile(abs.value, 0.25),xupper=quantile(abs.value, 0.75),
                   y = interaction(A_mod,  Y_mod,  M_mod), #interaction(A_fit,  M_fit, Y_fit), 
                 fill = Estimator #, shape = cluster_m
                 ), color = "black", size = 0.5, alpha = 1 , position = "dodge") +
  # geom_line(aes(y = value, x = sample_size, linetype = Estimator,
  #               group=interaction(Estimator,Nuisance_estimation),
  #                color = Nuisance_estimation) ) +
  scale_x_continuous("Relative Bias or RMSE" 
                     # ,breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)/10
                     ) +
  scale_y_discrete("fitted models for A, Y, M") +
  geom_vline(xintercept = 0.1, size = 0.2, linetype = "longdash") +
  # geom_vline(xintercept = 0.05, size = 0.2, linetype = "longdash") +
  # geom_hline(yintercept = 0, size = 0.4, linetype = "longdash") +
  # geom_hline(yintercept = -0.05, size = 0.2, linetype = "longdash") +
 # geom_hline(yintercept = -0.1, size = 0.2, linetype = "longdash") +
  # scale_x_continuous("Relative Bias", breaks = c(-c(9,5,3)/10, -0.1, 0.1, c(3,5,9)/10) ) +
  # scale_linetype_manual("Estimator", values = c("solid", "dashed")) +
  # scale_color_manual("Fit models", values = c("darkgreen", "blue")) +
  fig.theme + 
  theme(
    strip.background.x = element_rect(color = "lightgrey", fill = "lightgrey"),
    strip.background.y = element_blank(),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_blank(),
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15) # element_text(angle = 90, vjust = 0.5, hjust=1)
  )


fig.rbias

ggsave("Tables_Figs/fig_CLvsSL.pdf", height = 12.5, width = 12.5)





##  fixed vs RE ----
plotdf <- res_long1 %>% 
  filter(
    abs.performance %in% c("|Rel.Bias|", "RMSE"), 
    str_detect(effname, "Y\\(")
    # effname %in% c("Y(1,gm1(1),gm2(0))") # c( "Y(0,gmjo(0))" ), # (!true_model %in% c("other"))
    , J == 100, nj == 10 
    , cluster_a %in% c("FE","RE"), cluster_m %in% c("FE","RE"), cluster_y %in% c("FE","RE")
  )

fig.rbias <- ggplot(plotdf) +
  ggtitle( glue("{unique(plotdf$J)} clusters, cluster size {unique(plotdf$nj)}") ) + #, ICC = {unique(plotdf$icc)}, association of X and U = {unique(plotdf$x_z)};  estimand: {unique(plotdf$effname)}
  facet_grid(Estimator ~abs.performance , 
             #A_fit ~ Y_fit, #Y_true + Clus_true + Trt_true
             labeller = labeller(
               # cluster_a = label_both, cluster_y = label_both , 
               abs.performance = label_value, Estimator = label_both
             ), 
             scales = "free") +
  # geom_point(
  geom_boxplot(
    aes( x = abs.value, 
         # xmin = , xmax = max(abs.value), xlower=quantile(abs.value, 0.25),xupper=quantile(abs.value, 0.75),
         y = interaction(A_fit,  Y_fit, M_fit), #interaction(A_mod,  Y_mod,  M_mod), #
         fill = Estimator #, shape = cluster_m
    ), color = "black", size = 0.5, alpha = 1 , position = "dodge") +
  # geom_line(aes(y = value, x = sample_size, linetype = Estimator,
  #               group=interaction(Estimator,Nuisance_estimation),
  #                color = Nuisance_estimation) ) +
  scale_x_continuous("Relative Bias or RMSE" 
                     # ,breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)/10
  ) +
  scale_y_discrete("fitted models for A, Y, M") +
  geom_vline(xintercept = 0.1, size = 0.2, linetype = "longdash") +
  # geom_vline(xintercept = 0.05, size = 0.2, linetype = "longdash") +
  # geom_hline(yintercept = 0, size = 0.4, linetype = "longdash") +
  # geom_hline(yintercept = -0.05, size = 0.2, linetype = "longdash") +
  # geom_hline(yintercept = -0.1, size = 0.2, linetype = "longdash") +
  # scale_x_continuous("Relative Bias", breaks = c(-c(9,5,3)/10, -0.1, 0.1, c(3,5,9)/10) ) +
  # scale_linetype_manual("Estimator", values = c("solid", "dashed")) +
  # scale_color_manual("Fit models", values = c("darkgreen", "blue")) +
  fig.theme + 
  theme(
    strip.background.x = element_rect(color = "lightgrey", fill = "lightgrey"),
    strip.background.y = element_blank(),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_blank(),
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15) # element_text(angle = 90, vjust = 0.5, hjust=1)
  )

fig.rbias

ggsave("Tables_Figs/fig_FEvsRE.pdf", height = 12.5, width = 12.5)





## mse -----------------------

plotdf <- res_long1 %>% filter(
  performance %in% c("MSEsqrtN"), effect %in% c("meD", "reD"), (!true_model %in% c("other"))
) 


fig.mse <- ggplot(plotdf) +
  facet_grid(effect ~ true_model, #Y_true + Clus_true + Trt_true
             labeller = label_value, scales = "fixed") +
  geom_point(aes(y= value, # y = sqrtnMSE, # y=,
                 x = sample_size, shape = Estimator,
                 color = Nuisance_estimation), size = 2.5 ) +
  geom_line(aes(y = value, 
                x = sample_size, linetype = Estimator,
                group=interaction(Estimator,Nuisance_estimation),
                color = Nuisance_estimation) ) +
  scale_y_continuous("square root of n*MSE"  #"RMSE"
                     ) +
  # geom_hline(yintercept = 0.1, size = 0.4, linetype = "longdash") +
  # geom_hline(yintercept = 0, size = 0.4, linetype = "longdash") +
  # geom_hline(yintercept = -0.1, size = 0.4, linetype = "longdash") +
  # scale_x_continuous("Relative Bias", breaks = c(-c(9,5,3)/10, -0.1, 0.1, c(3,5,9)/10) ) +
  scale_linetype_manual("Estimator", values = c("solid", "dashed")) +
  # scale_color_manual("Fit models", values = c("darkgreen", "blue")) +
  fig.theme

fig.mse


ggsave("Tables_Figs/fig.sqrtnMSE.pdf", height = 10, width = 12)
# ggsave("Tables_Figs/fig.RMSE.pdf", height = 10, width = 12)


## coverage ------------
plotdf <- res_long1 %>% filter(
  performance %in% c("cover"), effect %in% c("meD", "reD"), (!true_model %in% c("other"))
)  


fig.cover <- ggplot(plotdf) +
  facet_grid(effect ~ true_model, #Y_true + Clus_true + Trt_true
             labeller = label_value, scales = "fixed") +
  geom_point(aes(y= value, # y = value, # y=abs.value,
                 x = sample_size, shape = Estimator,
                 color = Nuisance_estimation), size = 2.5 ) +
  geom_line(aes(y = value, 
                x = sample_size, linetype = Estimator,
                group=interaction(Estimator,Nuisance_estimation),
                color = Nuisance_estimation) ) +
  scale_y_continuous("Coverage Rate") +
  # geom_hline(yintercept = 0.1, size = 0.4, linetype = "longdash") +
  geom_hline(yintercept = 0.95, size = 0.4, linetype = "longdash") +
  # geom_hline(yintercept = -0.1, size = 0.4, linetype = "longdash") +
  # scale_x_continuous("Relative Bias", breaks = c(-c(9,5,3)/10, -0.1, 0.1, c(3,5,9)/10) ) +
  scale_linetype_manual("Estimator", values = c("solid", "dashed")) +
  # scale_color_manual("Fit models", values = c("darkgreen", "blue")) +
  fig.theme

fig.cover

ggsave("Tables_Figs/fig.cover.pdf", height = 10, width = 12)

# supplement 

## abs.bias -------

plotdf <- res_long1 %>% filter(
  performance %in% c("rbias"), effect %in% c("meD", "reD"), (!true_model %in% c("other"))
)  


fig.absrbias <- ggplot(plotdf) +
  facet_grid(effect ~ true_model, #Y_true + Clus_true + Trt_true
             labeller = label_value, scales = "fixed") +
  geom_point(aes(y = abs.value, 
                 x = sample_size, shape = Estimator,
                 color = Nuisance_estimation), size = 2.5 ) +
  geom_line(aes(y = abs.value, 
                x = sample_size, linetype = Estimator,
                group=interaction(Estimator,Nuisance_estimation),
                color = Nuisance_estimation) ) +
  scale_y_continuous("|Relative Bias|", breaks = c(0,0.5,1)/10 ) +
  geom_hline(yintercept = 0.1, size = 0.2, linetype = "longdash") +
  geom_hline(yintercept = 0, size = 0.4, linetype = "longdash") +
  # geom_hline(yintercept = -0.1, size = 0.4, linetype = "longdash") +
  # scale_x_continuous("Relative Bias", breaks = c(-c(9,5,3)/10, -0.1, 0.1, c(3,5,9)/10) ) +
  scale_linetype_manual("Estimator", values = c("solid", "dashed")) +
  # scale_color_manual("Fit models", values = c("darkgreen", "blue")) +
  fig.theme

fig.absrbias


ggsave("Tables_Figs/fig.absrbias.pdf", height = 10, width = 12)




save.image("/work/08878/xliu19/ls6/meMO/simulation/results.RData")


# NOT WORK ----------



# NOT USED currently --------

# Binomial
dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_Bqt_n300.RData" # updated quadratic.R means tt's true model is quadractic
jobconds <- which( (condition$quadratic.R== 1) & condition$Fit %in% c("mlr", "linear") & condition$n %in% c(300) )

dirpath <-"/work/08878/xliu19/ls6/meMO/simulation/RData/job_Bqt_n500.RData"
jobconds <- which( (condition$quadratic.R== 1) & condition$Fit %in% c("mlr", "linear") & condition$n %in% c(500) )


dirpath <-"/work/08878/xliu19/ls6/meMO/simulation/RData/job_Bqt_n1000.RData"
jobconds <- which( (condition$quadratic.R== 1) & condition$Fit %in% c("mlr", "linear") & condition$n %in% c(1000) )

dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_Ball_n300.RData"
jobconds <- which( (condition$quadratic.R + condition$quadratic.M + condition$quadratic.Y <= 3) & condition$Fit %in% c("mlr", "linear") & condition$n %in% c(300) )


dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_gauBq5000.RData"
jobconds <- which( (condition$quadratic.R + condition$quadratic.M + condition$quadratic.Y==3) & condition$Fit == "mlr" & condition$n %in% c(5000) )


dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_gauBq5000.RData"
jobconds <- which( (condition$quadratic.R + condition$quadratic.M + condition$quadratic.Y==3) & condition$Fit == "mlr" & condition$n %in% c(5000) ) 

dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_gauBq5000_Flin.RData"
jobconds <- which( (condition$quadratic.R + condition$quadratic.M + condition$quadratic.Y==3) & condition$Fit == "linear" & condition$n %in% c(5000) ) 


dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_gauBq1000_Fmlr.RData"
jobconds <- which( (condition$quadratic.R + condition$quadratic.M + condition$quadratic.Y==3) & condition$Fit == "mlr" & condition$n %in% c(1000) ) 


dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_gauBq1000_Fmlr_old.RData"
jobconds <- which( (condition$quadratic.R + condition$quadratic.M + condition$quadratic.Y==3) & condition$Fit == "mlr" & condition$n %in% c(1000) ) 

dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_gauBq1000_Flin_old.RData"
jobconds <- which( (condition$quadratic.R + condition$quadratic.M + condition$quadratic.Y==3) & condition$Fit == "linear" & condition$n %in% c(1000) ) 



dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_gauBqone1000_Flin_old.RData"
jobconds <- which( (condition$quadratic.R + condition$quadratic.M + condition$quadratic.Y < 3) & condition$Fit == "linear" & condition$n %in% c(1000) )


dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_gauBqone1000_Flin.RData"
jobconds <- which( (condition$quadratic.R + condition$quadratic.M + condition$quadratic.Y < 3) & condition$Fit == "linear" & condition$n %in% c(1000) )

dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_Ball_n1000_Fmlr.RData"
jobconds <- which( (condition$quadratic.R + condition$quadratic.M + condition$quadratic.Y <= 3) & condition$Fit == "mlr" & condition$n %in% c(1000) )



# 20 num_c
dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_Bqt_20c_n300_mlr.RData"
jobconds <- which( (condition$quadratic.R== 1) & condition$Fit %in% c("mlr") & condition$n %in% c(300) )

dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_Bnonqt_20c_n300_mlr.RData"
jobconds <- which( (condition$quadratic.R== 0) & condition$Fit %in% c("mlr") & condition$n %in% c(300) )



dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_Bqt_20c_n500_mlr2.RData"
jobconds <- which( (condition$quadratic.R== 1) & condition$Fit %in% c("mlr") & condition$n %in% c(500) )[2]

dirpath<-"/work/08878/xliu19/ls6/meMO/simulation/RData/old/job_Bqt_20c_n300_lin_allv.RData"
jobconds <- which( (condition$quadratic.R== 1) & condition$Fit %in% c("linear") & condition$n %in% c(300) )

dirpath<-"/work/08878/xliu19/ls6/meMO/simulation/RData/old/job_Bqt_20c_n300_mlr_allv.RData"
jobconds <- which( (condition$quadratic.R== 1) & condition$Fit %in% c("mlr") & condition$n %in% c(300) )


dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_Bqt_20c_n300_lin_safev.RData"
jobconds <- which( (condition$quadratic.R== 1) & condition$Fit %in% c("linear") & condition$n %in% c(300) )

#  "gaussian"
dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_Gqrqt_20c_n300_lin.RData"
jobconds <- which( ((condition$quadratic.R+condition$quadratic.tt) %in% c(1,2) ) & condition$Fit %in% c("linear") & condition$Yfamily=="gaussian" & condition$n %in% c(300) )




# "gaussian"
dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_gauN1000.RData"
jobconds <- which( (condition$quadratic.R + condition$quadratic.M + condition$quadratic.Y<=3) & condition$Fit == "mlr" & condition$n %in% c(1000) )


dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_gauN5000.RData"
jobconds <- which( (condition$quadratic.R + condition$quadratic.M + condition$quadratic.Y<=3) & condition$Fit == "mlr" & condition$n %in% c(5000) ) 

dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_gauNq5000.RData"
jobconds <- which( (condition$quadratic.R + condition$quadratic.M + condition$quadratic.Y==3) & condition$Fit == "mlr" & condition$n %in% c(5000) ) 


dirpath <- "/work/08878/xliu19/ls6/meMO/simulation/RData/job_gauNq1000.RData"
jobconds <- which( (condition$quadratic.R + condition$quadratic.M + condition$quadratic.Y==3) & condition$Fit == "mlr" & condition$n %in% c(1000) )










