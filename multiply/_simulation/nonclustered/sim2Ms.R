library(glue)
library(tidyverse)
library(mvtnorm)
library(SuperLearner)
# library(devtools)
# load_all("multiplymed")

source("_simulation/clustered/gen_clustered.R")

source("_simulation/SL.lightgbm.R")
# source("_simulation/SL.glm.saturated.R")
# source("_simulation/SL.glmnet3.R")


OneData <- function(
        ) {
    dat <- gendata(N = 1000)
    Wnames <- names(dat)[startsWith(names(dat), "W")] 
    Aname <- "A"
    Yname <- "Y"
    Znames <- names(dat)[startsWith(names(dat), "Z")]
    M1names <- names(dat)[startsWith(names(dat), "M1")]
    M2names <- names(dat)[startsWith(names(dat), "M2")]
    
    mediation(dat, "A", "W1", 
              z, m, "Y", 
              family = "binomial", 
              folds = folds, partial_tmle = tmle,
              learners_g = "SL.mean", 
              learners_c = learners,
              learners_b = learners, 
              learners_e = learners, 
              learners_hz = learners, 
              learners_u = learners, 
              learners_ubar = learners, 
              learners_v = learners, 
              learners_vbar = learners)
} 

res <- mutate(res, 
              n = case_when(
                  n == 1 ~ 500, 
                  n == 2 ~ 1000, 
                  n == 3 ~ 5000, 
                  n == 4 ~ 1e4
              ))

saveRDS(res, glue("_research/data/sim_not_transported_{tmle}_{dgp}_{id}.rds"))
