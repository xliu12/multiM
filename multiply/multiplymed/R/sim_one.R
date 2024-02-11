library(mvtnorm)
library(lme4)
# library(cubature)
library(parallel)
library(tidyverse)
library(data.table)
library(fastDummies)
library(glue)
library(SuperLearner)
library(origami)


devtools::load_all("multiplymed")
# source("./multiplymed/R/fit_pred.R")
# source("./multiplymed/R/y.R")
# source("./multiplymed/R/am.R")
# source("./multiplymed/R/twoMcl.R")

# source("multiplymed/R/fun_gendata.R")


GenData <- function(
        seedone = 123,
        num_clust = 100,
        clust_size = 30,
        quadratic.A = FALSE,
        quadratic.M = FALSE,
        quadratic.Y = FALSE,
        num_m = 2,
        num_x = 3,
        Yfamily = "gaussian"
) {
   
    N <- num_clust * clust_size
    
    data <- data.frame(
        id = 1:N, 
        school = rep(1:num_clust, each = clust_size)
    )
    
    # X individual-level confounders ------
    iccx <- 0.2
    xb <- mvtnorm::rmvnorm(n = num_clust, 
                           mean = rep(0, num_x), 
                           sigma = diag(iccx, nrow = num_x))[rep(1:num_clust, each = clust_size), ]
    xe <- mvtnorm::rmvnorm(n = N, 
                           mean = rep(0, num_x), 
                           sigma = diag(1-iccx, nrow = num_x))
    x <- xb + xe
    data$X <- x
    
    # Z (unobserved) cluster-level confounders -----------------
    z <- rep(rnorm(num_clust, sd = 1), each = clust_size)
    data$Z <- z
    
    # Treatment ---------
    gen_a <- list(icca = 0.2, a_x = sqrt(0.15 * 1 / num_x), a_z = sqrt(0.15 / 1))
    ab <- rep(rnorm(num_clust, sd = sqrt(gen_a[["icca"]])), each = clust_size) 
    if (quadratic.A == FALSE) {
        a_given <- ab +  gen_a[["a_x"]] * rowSums(data$X) + gen_a[["a_z"]] * data$Z
    }
    if (quadratic.A == TRUE) {
        Xquad <- (data$X ^ 2 - 1) / sqrt(2) # mean(rnorm(.)^2) = 1; var(rnorm(.)^2) = 2
        a_given <- ab +  gen_a[["a_x"]] * rowSums(Xquad) + gen_a[["a_z"]] * data$Z
    }
    data$A <- rbinom(N, 1, pa(1, a_given, gen_a))
    
    # Mediators ------------------
    gen_m <- list(iccm1 = 0.2, iccm2 = 0.2, 
                  m1_on_a = 0.3, m1_on_x = sqrt(0.15 / num_x), m1_on_z = sqrt(0.15), 
                  m2_on_a = 0.2, m2_on_m1 = 0.2, m2_on_am1 = 0.3,
                  m2_on_x = sqrt(0.15 / num_x),
                  m2_on_z = sqrt(0.15) 
    )
    
    m1b <- rep(rnorm(num_clust, sd = sqrt(gen_m[["iccm1"]])), each = clust_size)
    m2b <- rep(rnorm(num_clust, sd = sqrt(gen_m[["iccm2"]])), each = clust_size)

    if (quadratic.M == TRUE) {
        Xquad <- (data$X ^ 2 - 1) / sqrt(2) 
        m1_given <- m1b + gen_m[["m1_on_x"]] * rowSums(Xquad) + 
            gen_m[["m1_on_z"]] * data$Z
        data$M1 <- rbinom(N, 1, prob = pm1(1, data$A, m1_given, gen_m) )
        
        m2_given <- m2b + gen_m[["m2_on_x"]] * rowSums(Xquad) +
            gen_m[["m2_on_z"]] * data$Z
        data$M2 <- rbinom(N, 1, prob = pm2(1, data$M1, data$A, m2_given, gen_m) )
        
    }
    if (quadratic.M == FALSE) {
        m1_given <- m1b + gen_m[["m1_on_x"]] * rowSums(data$X) + 
            gen_m[["m1_on_z"]] * data$Z
        # pm1(1, data$A, m1_given, gen_m) %>% quantile(.,c(0.1, 0.9))
        data$M1 <- rbinom(N, 1, prob = pm1(1, data$A, m1_given, gen_m) )
        
        m2_given <- m2b + gen_m[["m2_on_x"]] * rowSums(data$X) +
            gen_m[["m2_on_z"]] * data$Z
        data$M2 <- rbinom(N, 1, prob = pm2(1, data$M1, data$A, m2_given, gen_m) )
    }
    
    # Outcome -------------------
    
    gen_y <- list(iccy = 0.2, yintercept = 1,
                  y_on_a = 0.3, y_on_m1 = 0.3, y_on_m2 = 0.3,
                  y_on_am1 = 0, y_on_am2 = 0, y_on_m1m2 = 0.6, y_on_am1m2 = 0,
                  y_on_x = sqrt(0.35 / num_x), y_on_z = sqrt(0.15)
                  )
    yb <- rnorm(num_clust, sd = sqrt(gen_y[["iccy"]]))[rep(1:num_clust, each = clust_size)] 
    
    if (quadratic.Y == TRUE) {
        Xquad <- (data$X ^ 2 - 1) / sqrt(2)
        y_given <- yb + gen_y[["yintercept"]] + gen_y[["y_on_x"]] * rowSums(Xquad) +
            gen_y[["y_on_z"]] * data$Z
    }
    if (quadratic.Y == FALSE) {
        y_given <- yb + gen_y[["yintercept"]] + gen_y[["y_on_x"]] * rowSums(data$X) +
            gen_y[["y_on_z"]] * data$Z
    }
    
    if (Yfamily == "gaussian") {
        condmy <- my(m2 = data$M2, m1 = data$M1, data$A, given = y_given, gen_y, binary = FALSE)
        data$Y <- condmy + rnorm(N, sd = sqrt(1 - gen_y[["iccy"]])) 
    }
    if (Yfamily == "binomial") {
        condmy <- my(m2 = data$M2, m1 = data$M1, data$A, given = y_given, gen_y, binary = TRUE)
        data$Y <- rbinom(N, 1, condmy)
    }
    
    datobs <- do.call(data.frame, data)
    
    
    trueVals <- function() {
        a_vals12 <- rbind(
            # IIE_M1, fixing A=1 and p(M2(0) | c)
            expand.grid(a0 = c(1), a1 = c(0, 1), a2 = c(0)),
            # IIE_M2, fixing A=1 and p(M1(1) | c)
            expand.grid(a0 = c(1), a1 = c(1), a2 = c(0, 1))
        )
        a_vals12$ajo <- NA
        a_vals <- rbind(a_vals12, 
                        data.frame(a0 = c(1, 1, 0), 
                                   a1 = NA, a2 = NA,
                                   ajo = c(0, 1, 0)))
        truevals <- list()
        for (j in 1:nrow(a_vals)) {
            a0 <- a_vals$a0[j]
            a1 <- a_vals$a1[j]
            a2 <- a_vals$a2[j]
            ajo <- a_vals$ajo[j]
            
            vals <- expand.grid(m1 = c(0, 1), m2 = c(0, 1))
            if (!is.na(a1)) {
                # integrate over the marginal mediator distributions, M1 under a1 and M2 under a2 
                intm1m2 <- sapply(1:nrow(vals), function(i=1) {
                    m1 <- vals$m1[i]
                    m2 <- vals$m2[i]
                    p_marginal <- pm1(m1, a1, m1_given, gen_m) * 
                        pm2a(m2, a2, m2_given, gen_m) 
                    
                    p_marginal * my(m2 = vals$m2[i], m1 = vals$m1[i], a = a0, 
                           y_given, gen_y, binary = (Yfamily=="binomial")) 
                    
                }, simplify = TRUE)
                int_m1m2 <- rowSums(intm1m2)
                truevals[[glue("Y({a0},gm1({a1}),gm2({a2}))")]] <- mean(int_m1m2)
            }
            
            if (!is.na(ajo)) {
                # integrate over the joint mediator distribution, (M1, M2) under ajo
                intm1m2 <- sapply(1:nrow(vals), function(i=1) {
                    m1 <- vals$m1[i]
                    m2 <- vals$m2[i]
                    p_joint <- pm1(m1, ajo, m1_given, gen_m) * 
                        pm2(m2, m1, ajo, m2_given, gen_m) 
                    
                    p_joint * my(m2 = vals$m2[i], m1 = vals$m1[i], a = a0, 
                       y_given, gen_y, binary = (Yfamily=="binomial")) 
                    
                }, simplify = TRUE)
                int_m1m2 <- rowSums(intm1m2)
                truevals[[glue("Y({a0},gmjo({ajo}))")]] <- mean(int_m1m2)
            }
        }
        
        return(truevals)
    }
    
    truevals <- trueVals()
    
    out <- mget(ls(envir = environment()))
    
    return(out)
    
}






# Simulation----------------------------

condition_all <- data.frame(expand.grid(
    num_clust = c(100),
    clust_size = c(100),
    # generating data
    quadratic.A = c(F, T),
    quadratic.M = c(F, T),
    quadratic.Y = c(F, T),
    Yfamily = c("gaussian"), # c("gaussian"),
    cluster_a = "FE", # "RE", # "noncluster", # 
    cluster_m =  "FE", # "RE", #  "noncluster.mlr", # 
    cluster_y =  "FE", # "noncluster.mlr", ## "FE.mlr", # 
    # interaction_fitm2 =  c("AM.M1"), # NULL, # 
    # interaction_fity = c("AM.M1", "AM.M2", "Mint", "AMint"), # NULL, # 
    Fit = c("glm") # , "mlr"
)) 

condition <- condition_all[(condition_all$quadratic.A + condition_all$quadratic.M + condition_all$quadratic.Y) %in% c(1, 3, 0), ]

condition <- condition %>% 
    mutate(
        cluster_opt_a = paste(cluster_a, Fit, sep="."),
        cluster_opt_m = paste(cluster_m, Fit, sep="."),
        cluster_opt_y = paste(cluster_y, Fit, sep=".")
        )


set.seed(12)
datseeds <- sample(1:1e6, 1000)

iseed <-1
cond <- 1 



OneData <- function(iseed = 1, cond = 1){
    Yfamily <- as.character(condition$Yfamily[cond])
    
    gen_data <- GenData(
        seedone = datseeds[iseed], 
        num_clust = condition$num_clust[cond],
        clust_size = condition$clust_size[cond],
        quadratic.A = condition$quadratic.A[cond],
        quadratic.M = condition$quadratic.M[cond],
        quadratic.Y = condition$quadratic.Y[cond],
        Yfamily = Yfamily,
        num_m = 2,
        num_x = 3
    )
    
    data <- gen_data$datobs
    Sname <- "school"
    Xnames <- colnames(data)[grep("^X", colnames(data))]
    Wnames <- colnames(data)[grep("^W", colnames(data))]
    
    Aname <- "A"
    Mnames <- colnames(data)[grep("^M", colnames(data))]
    Yname <- "Y"
    
    # true values -----------------------------
    gen_largeN <- GenData(
        seedone = datseeds[iseed], 
        num_clust = 600,
        clust_size = 500,
        quadratic.A = condition$quadratic.A[cond],
        quadratic.M = condition$quadratic.M[cond],
        quadratic.Y = condition$quadratic.Y[cond],
        Yfamily = Yfamily,
        num_m = 2,
        num_x = 3
    )
    true_values <- unlist(gen_largeN$truevals)
    
    true_effs <- effect(true_values)
    
    # gen_data$gen_m$m2_on_am1
    # gen_data$gen_y$y_on_m1m2
    # gen_data$quadratic.A; gen_data$quadratic.M; gen_data$quadratic.Y
    
    # estimators -----
    Fit <- as.character(condition$Fit[cond])
    if (Fit == "glm") { 
        learners_a <- learners_m <- learners_y <- c("SL.glm") 
    }
    if (Fit == "mlr") { 
        learners_a <- learners_m <- learners_y <- c("SL.glm", "SL.gam", "SL.lightgbm") 
    }
    cluster_opt_a <- condition$cluster_opt_a[cond]
    cluster_opt_m <- condition$cluster_opt_m[cond]
    cluster_opt_y <- condition$cluster_opt_y[cond]
    
    out <- clustered(data, 
              Sname,
              Wnames = NULL, Xnames,
              Aname,
              Mnames, Mfamily = "binomial",
              Yname, Yfamily,
              cluster_opt_a = cluster_opt_a, 
              # "FE.mlr", # "RE", # "noncluster.mlr", # 
              cluster_opt_m = cluster_opt_m, # "FE.mlr", # "RE", #  "noncluster.mlr", # 
              cluster_opt_y = cluster_opt_y, # "FE.mlr", # "noncluster.mlr", ## "FE.mlr", # 
              interaction_fitm2 =  c("AM.M1"), # NULL, # 
              interaction_fity = c("AM.M1", "AM.M2", "Mint", "AMint"), # NULL, # 
              learners_a = learners_a, # c("SL.glm", "SL.gam"),
              learners_m = learners_m, #c("SL.glm", "SL.gam"),
              learners_y = learners_y, #c("SL.glm", "SL.gam")
              )
    
    # unlist(out$thetas) / true_values - 1
    # unlist(out$regs) / true_values - 1
    # unlist(out$rmpw) / true_values - 1

    estimates <- cbind(
        multi = out$thetas_effs,
        se_multi = out$se_effs,
        ci1_multi = out$ci_effs[1, ],
        ci2_multi = out$ci_effs[2, ],
        reg = out$regs_effs,
        rmpw = out$rmpw_effs)
    
    out <- data.frame(condition[cond, ],
                      effname = names(true_effs), 
                      true_val = true_effs,
                      estimates)
    
    return(out)
}


# Parallel setup -----
num_reps <- 20 
mc_cores <- 20
seedseq <- seq(from=1, to=1000, by=num_reps)

jobconds <- c(1:nrow(condition))

# zip file to upload -----
zip("multiplymed.zip", files = dir("multiplymed"))
