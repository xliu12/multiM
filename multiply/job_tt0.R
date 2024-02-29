library(parallel)

# setwd("/work/08878/xliu19/ls6/multiply")

devtools::load_all("multiplymed")
# source("multiplymed/R/fit_pred.R")
# source("multiplymed/R/y.R")
# source("multiplymed/R/am.R")
# source("multiplymed/R/twoMcl.R")
source("simulation/sim_one.R")
source("simulation/fun_gendata.R")

# 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

# Simulation conditions ----


set.seed(12)

# condition <- condition %>% 
#   # filter(quadratic.A + quadratic.M + quadratic.Y == 0) %>%
#   filter(cluster_y == "noncluster")

jobconds <- c(1:nrow(condition))

reslist <- list()

jobseeds <- 1:1000 

nj <- paste(unique(condition$clust_size), collapse = "-")
icc <- paste(unique(condition$icc)*10, collapse = "-")
J <- paste(unique(condition$num_clust), collapse = "-")
learner <- unique(condition$Fit)
Yfamily <- ifelse(unique(condition$Yfamily) == "binomial", "bin", "linear")
xz <- paste(unique(condition$x_z)*10, collapse = "-")
rname <- glue("RData/xz{xz}_{Yfamily}_Fit_{learner}_icc{icc}_J{J}_n{nj}_rep{jobseeds[1]}_{tail(jobseeds, 1)}_addmethods.RData")
rname <- as.character(rname)


for(k in jobconds) {
  res <- mclapply(jobseeds, 
                  OneData, cond = k, 
                  mc.preschedule = TRUE, mc.cores = 100) 
  reslist[[k]] <- res
  save(
    reslist, file = as.character(rname) # paste0("RData/linear_Fit_glm_nj50", ".RData")
  )
}


