library(mvtnorm)
library(lme4)
# library(cubature)
library(parallel)
library(tidyverse)
library(data.table)
library(fastDummies)
library(glue)
library(SuperLearner)
library(ranger)
library(lightgbm)
library(origami)

devtools::load_all("multiplymed")

load("_example/example_data.RData")

# Estimation ----------------------

colnames(datobs)
Sname <- "school" 
Xnames <- colnames(datobs) %>% grep("^X_", ., value = TRUE)
Wnames <- c("W_BYURBAN.x_2", "W_BYACCLIM", "W_BY10FLP") 

Aname <- "A" 
Mnames <- c("M1", "M2") 
Yname <- "Y" 



methds <- data.frame(
  cluster_a = c("FE"), 
  cluster_m =  c("FE"), 
  cluster_y =  c("FE"), 
  Fit = c("mlr") 
)


set.seed(12)
datseeds <- sample(1:1e6, 1000)

iseed <-0
cond <- "M1M2"


OneData <- function(iseed = 1, cond = "M1M2") {
  if (iseed != 0) {
    set.seed(datseeds[iseed])
    J <- length(unique(datobs$school))
    boot_Sj <- sample(unique(datobs$school), J, replace = TRUE)
    dat_boot <- datobs %>% 
      filter(school %in% boot_Sj)
  }
  if (iseed ==0) {
    dat_boot <- datobs
  }
  
  
  if (cond == "M1M2") {
    Mnames <- c("M1", "M2")
  }
  if (cond == "M2M1") {
    Mnames <- c("M2", "M1")
  }
  
  if(length(unique(datobs$Y)) == 2) {Yfamily <- "binomial"}
  if(length(unique(datobs$Y)) > 2) {Yfamily <- "gaussian"}
  data <- dat_boot
  
  one.jj <- function(jj = 1) {
    Fit <- as.character(methds$Fit[jj])
    if (Fit == "glm") {
      learners_a <- learners_m <- learners_y <- c("SL.glm")
      num_folds <- 1
    }
    if (Fit == "mlr") {
      learners_a <- learners_m <- learners_y <- c("SL.ranger", "SL.glm")
      Xnames <- c(Xnames, Wnames)
      num_folds <- 1
    }
    
    cluster_opt_a <- methds$cluster_opt_a[jj]
    cluster_opt_m <- methds$cluster_opt_m[jj]
    cluster_opt_y <- methds$cluster_opt_y[jj]
    
    out <- clustered(data,
                     Sname,
                     Wnames, # = NULL, 
                     Xnames,
                     Aname,
                     Mnames, Mfamily = "binomial",
                     Yname, Yfamily,
                     cluster_opt_a = cluster_opt_a,
                     cluster_opt_m = cluster_opt_m, 
                     cluster_opt_y = cluster_opt_y, 
                     interaction_fitm2 =  c("AM.M1"), 
                     interaction_fity = c("AM.M1", "AM.M2", "Mint", "AMint"), 
                     num_folds = num_folds,
                     learners_a = learners_a, 
                     learners_m = learners_m, 
                     learners_y = learners_y, 
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


# run -----------------------
jobseeds <- 0:1000

rname <- "example_results.RData"
rname <- as.character(rname)

res <- mclapply(jobseeds, 
                OneData, cond = "M1M2", 
                mc.preschedule = TRUE, mc.cores = 4) 

save(
  res, file = as.character(rname)  
)




