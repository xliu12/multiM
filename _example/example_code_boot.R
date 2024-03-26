library(mvtnorm)
library(lme4)
library(parallel)
library(tidyverse)
library(data.table)
library(fastDummies)
library(glue)
library(SuperLearner)
library(ranger)
library(lightgbm)
library(origami)

devtools::load_all("MediationClustered")

load("_example/example_data.RData")
dat <- data_example

# Estimation ----------------------


Sname <- "school" 
Xnames <- colnames(dat) %>% grep("^X_", ., value = TRUE)
Aname <- "A" 
Mnames <- c("M1", "M2") 
Yname <- "Y" 

methds <- data.frame(
  cluster_a = c("FE"), 
  cluster_m =  c("FE"), 
  cluster_y =  c("FE"), 
  Fit = c("mlr") 
) %>%
  mutate(
    cluster_opt_a = paste(cluster_a, Fit, sep="."),
    cluster_opt_m = paste(cluster_m, Fit, sep="."),
    cluster_opt_y = paste(cluster_y, Fit, sep=".")
  ) 


set.seed(12)
datseeds <- sample(1:1e6, 1000)

iseed <-0


OneData <- function(iseed = 1) {
  if (iseed != 0) {
    set.seed(datseeds[iseed])
    J <- length(unique(dat$school))
    boot_Sj <- sample(unique(dat$school), J, replace = TRUE)
    dat_boot <- dat %>% 
      filter(school %in% boot_Sj)
  }
  if (iseed ==0) {
    dat_boot <- dat
  }
  
  
  if(length(unique(dat$Y)) == 2) {Yfamily <- "binomial"}
  if(length(unique(dat$Y)) > 2) {Yfamily <- "gaussian"}
  data <- dat_boot
  
  one.jj <- function(jj = 1) {
    Fit <- as.character(methds$Fit[jj])
    if (Fit == "glm") {
      learners_a <- learners_m <- learners_y <- c("SL.glm")
    }
    if (Fit == "mlr") {
      learners_a <- learners_m <- learners_y <- c("SL.ranger", "SL.earth")
    }
    
    cluster_opt_a <- methds$cluster_opt_a[jj]
    cluster_opt_m <- methds$cluster_opt_m[jj]
    cluster_opt_y <- methds$cluster_opt_y[jj]
    
    out <- clustered(data,
                     Sname = Sname,
                     Wnames = NULL, 
                     Xnames = Xnames,
                     Aname = Aname,
                     Mnames = Mnames, 
                     Yname = Yname, 
                     Yfamily = Yfamily,
                     cluster_opt_a = cluster_opt_a,
                     cluster_opt_m = cluster_opt_m, 
                     cluster_opt_y = cluster_opt_y, 
                     interaction_fitm2 =  c("AM.M1"), 
                     interaction_fity = c("AM.M1", "AM.M2", "Mint", "AMint"), 
                     learners_a = learners_a, 
                     learners_m = learners_m, 
                     learners_y = learners_y
    )
    
    
    estimates <- cbind(
      est_multi = out$thetas_effs
      )
    
    res <- data.frame(methds[jj, ],
                      effname = names(out$thetas_effs),
                      estimates, row.names = NULL)
    
    return(res)
  }
  
  oneboot <- lapply(1:nrow(methds), one.jj)
  one_boot <- do.call(rbind, oneboot)
  
  return(one_boot)
}


# run 
jobseeds <- 1:1000

rname <- "example_results.RData"
rname <- as.character(rname)

res <- mclapply(jobseeds, 
                OneData, 
                mc.preschedule = TRUE, mc.cores = 4) 

save(
  res, file = as.character(rname)  
)




