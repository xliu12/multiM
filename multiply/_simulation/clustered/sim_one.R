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

source("./_simulation/clustered/fun_gendata.R")

source("./multiplymed/R/fit_pred.R")
source("./multiplymed/R/y.R")
source("./multiplymed/R/am.R")
source("./multiplymed/R/twoMcl.R")


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

condition1 <- data.frame(expand.grid(
    num_clust = c(100),
    clust_size = c(100),
    # generating data
    quadratic.A = c(F, T),
    quadratic.M = c(F, T),
    quadratic.Y = c(F, T),
    Yfamily = c("gaussian"), # c("gaussian"),
    Fit =  c("linear", "mlr")
)) 
adding1 <- data.frame(expand.grid(
    num_clust = c(100),
    clust_size = c(100),
    # generating data
    quadratic.A = c(F, T),
    quadratic.M = c(F, T),
    quadratic.Y = c(F, T),
    Yfamily = c("binomial"),  
    Fit =  c("linear", "mlr")
))
condition_all <- data.frame(do.call(rbind, list(condition1, adding1)))
condition <- condition_all
#[(condition_all$quadratic.A + condition_all$quadratic.M + condition_all$quadratic.Y) %in% c(1, 3, 0), ]




set.seed(12)
datseeds <- sample(1:1e6, 1000)

iseed <-1
cond <- 16 # 100        100        TRUE        TRUE        TRUE gaussian    



OneData <- function(iseed = 1, cond = 1){
    Yfamily <- "gaussian" #  as.character(condition$Yfamily[cond])
    
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
        num_clust = 1000,
        clust_size = 2e3,
        # quadratic.A = condition$quadratic.A[cond],
        # quadratic.M = condition$quadratic.M[cond],
        # quadratic.Y = condition$quadratic.Y[cond],
        Yfamily = Yfamily,
        num_m = 2,
        num_x = 3
    )
    true_values <- unlist(gen_largeN$truevals)
    
    effect(true_values)
    
    gen_data$gen_m$m2_on_am1
    gen_data$gen_y$y_on_m1m2
    gen_data$quadratic.A; gen_data$quadratic.M; gen_data$quadratic.Y
    
    # estimators -----
    Fit <- as.character(condition$Fit[cond])
    
    out <- clustered(data, 
              Sname,
              Wnames = NULL, Xnames,
              Aname,
              Mnames, Mfamily = "binomial",
              Yname, Yfamily,
              cluster_opt_a = "FE.mlr", # "RE", # "noncluster.mlr", # 
              cluster_opt_m =  "FE.mlr", # "RE", #  "noncluster.mlr", # 
              cluster_opt_y =  "FE.mlr", # "noncluster.mlr", ## "FE.mlr", # 
              interaction_fitm2 =  c("AM.M1"), # NULL, # 
              interaction_fity = c("AM.M1", "AM.M2", "Mint", "AMint"), # NULL, # 
              learners_a = c("SL.glm", "SL.gam"),
              learners_m = c("SL.glm", "SL.gam"),
              learners_y = c("SL.glm", "SL.gam")
              )
    
    unlist(out$thetas) / true_values - 1
    unlist(out$regs) / true_values - 1
    unlist(out$rmpw) / true_values - 1

    
    # crossfit_out$estimates / true_effects - 1
    # tmle_out$tml_estimates /true_effects-1
    
    out <- list(
        results = data.frame(
            condition[cond, ],
            effect = c("Yt1r1", "Yt1r0", "Yt1r1.Mt1r0", 
                       "Yt0r1", "Yt0r0", "Yt0r1.Mt0r0", 
                       "reD",  
                       "meD", 
                       "toD"),
            # true
            true_val = t(true_effects),
            # tml
            tml_estimate = tmle_out$tml_estimates,
            tml_interval = t(tmle_out$tml_intervals),
            # onestep
            onestep_estimate = crossfit_out$estimates,
            onestep_interval = t(crossfit_out$z_intervals),
            stderr = crossfit_out$stderrors
            # weighting-based
            , weighting_estimate = crossfit_out$weighting_estimates
        ) #, 
        # crossfit_out = crossfit_out, 
        # tmle_out = tmle_out
    )
    
    return(out)
}


# Parallel setup -----
num_reps <- 20 
mc_cores <- 20
seedseq <- seq(from=1, to=1000, by=num_reps)

jobconds <- c(1:nrow(condition))

# w_n <- function(n = 200, num_w = 3) { # w: baseline covariates
#     rmvnorm(n, mean = rep(0, num_w))
# }
# g <- function(a, w, randomized = FALSE) {
#     pscore <- pnorm( rowSums(sqrt(0.13/ncol(w)) * w) )
#     if (randomized) {pscore <- .5} 
#     a * pscore + (1 - a) * (1 - pscore)
# }
# 
# pm <- function(m1, m2, z, a, w, binary = TRUE) {
#     U <- rnorm(nrow(w)) # common causes of mediators
#     Em1 <- 0.1*z + 0.1*a + 0.2*U + rowSums(sqrt(0.13/ncol(w)) * w)
#     Em2 <- 0.3*z + 0.3*a + 0.2*U + rowSums(sqrt(0.13/ncol(w)) * w)
#     
#     if (!binary) {
#         em <- rmvnorm(nrow(w), mean = c(0,0), sigma = diag((1-0.2^2-0.13), nrow = 2))
#         out <- cbind(Em1, Em2) + em
#     }
#     if (binary) { 
#         prob1 <- pnorm(Em1)
#         prob2 <- pnorm(Em2)
#         
#         out <- cbind(m1 * prob1 + (1 - m1) * (1 - prob1),
#                      m2 * prob2 + (1 - m2) * (1 - prob2) )
#     }
#     out
# }
# 
# pm1aw <- function(m1, a, w) {
#     pz(1, a, w) * (
#         apply(pm(m1, m2=1, 1, a, w), 1, prod) + 
#             apply(pm(m1, m2=0, 1, a, w), 1, prod) )  +
#         pz(0, a, w) * (
#             apply(pm(m1, m2=1, 0, a, w), 1, prod) + 
#                 apply(pm(m1, m2=0, 0, a, w), 1, prod) ) 
# }
# 
# pm2aw <- function(m2, a, w) {
#     pz(1, a, w) * (
#         apply(pm(m1=1, m2, 1, a, w), 1, prod) + 
#             apply(pm(m1=0, m2, 1, a, w), 1, prod) )  +
#         pz(0, a, w) * (
#             apply(pm(m1=1, m2, 0, a, w), 1, prod) + 
#                 apply(pm(m1=0, m2, 0, a, w), 1, prod) ) 
# }
# 
# pm1m2aw <- function(m1, m2, a1, a2, w) {
#     pm1aw(m1, a1, w) * pm2aw(m2, a2, w)
# }
# 
# pmaw <- function(m1, m2, a, w) {
#     pz(1, a, w) * apply(pm(m1, m2, 1, a, w), 1, prod) +
#         pz(0, a, w) * apply(pm(m1, m2, 0, a, w), 1, prod)
# }
# 
# 
# m_ratio <- function(m1, m2, z, a, w) {
#     pm1aw(m1, a, w) * pm2aw(m2, a, w) / apply(pm(m1, m2, z, a, w), 1, prod)
# }
# 
# 
# 
# my <- function(m1, m2, z, a, w, binary = TRUE) {
#     mzaw <- 0.3*m1 + 0.3*m2 + 0.3*z + 0.3*a + rowSums(sqrt(0.13/ncol(w)) * w) # no interaction
#     if (binary) {
#         out <- pnorm(mzaw)
#     }
#     if (!binary) {
#         out <- mzaw + rnorm(nrow(w), sd = sqrt(1-0.13))
#     }
#     out
# }
# 
# u <- function(z, w, aprime, astar) {
#     my(1, 1, z, aprime, w) * pmaw(1, 1, astar, w) +
#         my(0, 1, z, aprime, w) * pmaw(0, 1, astar, w) +
#         my(0, 0, z, aprime, w) * pmaw(0, 0, astar, w) +
#         my(1, 0, z, aprime, w) * pmaw(1, 0, astar, w)
# }
# 
# intu <- function(w, aprime, astar) {
#     u(1, w, aprime, astar) * pz(1, aprime, w) +
#         u(0, w, aprime, astar) * pz(0, aprime, w) 
# }
# 
# intv <- function(m1, m2, w, aprime) {
#     my(m1, m2, 1, aprime, w) * pz(1, aprime, w) +
#         my(m1, m2, 0, aprime, w) * pz(0, aprime, w) 
# }
# 
# 
# 
# gendata <- function(N = 1000) {
#     # w1 <- rbinom(N, 1, .4)
#     # w2 <- rbinom(N, 1, .4)
#     w <- w_n(N, 3)
#     
#     a <- rbinom(N, 1, g(1, w))
#     z <- rbinom(N, 1, pz(1, a, w))
#     pm_jo <- pm(1, 1, z, a, w)
#     m1 <- rbinom(N, 1, pm_jo[, 1])
#     m2 <- rbinom(N, 1, pm_jo[, 2])
#     y <- rbinom(N, 1, my(m1, m2, z, a, w))
#     data.frame(W1 = w, A = a, Z1 = z, M1 = m1, M2 = m2, Y = y)
# }
# 
# 
# 
# 
# truth <- function(Nsuper = 1e4) {
#     w <- w_n(Nsuper, 3)
#     
#     
#     true_thetas <- rep(NA, 6)
#     names(true_thetas) <- c("Y(0,M(t_jo=0))", "Y(1,M(t_jo=0))", "Y(1,M(t_jo=1))", 
#                             "Y(1,M1(1),M2(1))", "Y(1,M1(0),M2(1))", "Y(1,M1(0),M2(0))")
#     
#     for (param in list(c(t0 = 0, t_jo = 0, t_1 = NA, t_2 = NA), 
#                        c(t0 = 1, t_jo = 0, t_1 = NA, t_2 = NA), 
#                        c(t0 = 1, t_jo = 1, t_1 = NA, t_2 = NA),
#                        c(t0 = 1, t_jo = NA, t_1 = 1, t_2 = 1),
#                        c(t0 = 1, t_jo = NA, t_1 = 0, t_2 = 1),
#                        c(t0 = 1, t_jo = NA, t_1 = 0, t_2 = 0)
#     )) {
#         t0 <- param["t0"]
#         if (!is.na(param["t_jo"])) {
#             t_jo <- param["t_jo"]
#             
#             true_v <- intv(1, 1, w, t0) * pmaw(1, 1, t_jo, w) +
#                 intv(1, 0, w, t0) * pmaw(1, 0, t_jo, w) +
#                 intv(0, 1, w, t0) * pmaw(0, 1, t_jo, w) +
#                 intv(0, 0, w, t0) * pmaw(0, 0, t_jo, w)
#             
#             vname <- paste0("Y(", t0 ,",M(t_jo=", t_jo, "))")
#             true_thetas[vname] <- mean(true_v)
#         }
#         if (is.na(param[["t_jo"]])) {
#             t_1 <- param["t_1"]
#             t_2 <- param["t_2"]
#             
#             true_v <- intv(1, 1, w, t0) * pm1m2aw(1, 1, t_1, t_2, w) +
#                 intv(1, 0, w, t0) * pm1m2aw(1, 0, t_1, t_2, w) +
#                 intv(0, 1, w, t0) * pm1m2aw(0, 1, t_1, t_2, w) +
#                 intv(0, 0, w, t0) * pm1m2aw(0, 0, t_1, t_2, w)
#             
#             vname <- paste0("Y(", t0 ,",M1(", t_1, "),", "M2(", t_2, ")", ")")
#             # true_thetas[vname] <- weighted.mean(true_v, prob_w)
#             true_thetas[vname] <- mean(true_v)
#         }
#     }
#     
#     true_effs <- effect(true_thetas)
#     
#     true_effs
# }
# 
# 
# # trueval <- truth()
# 
# If <- function(dat) {
#     w <- dat[, Wnames, drop = F]
#     true_Ifs <- rep(NA, 6)
#     names(true_Ifs) <- c("Y(0,M(t_jo=0))", "Y(1,M(t_jo=0))", "Y(1,M(t_jo=1))", 
#                          "Y(1,M1(1),M2(1))", "Y(1,M1(0),M2(1))", "Y(1,M1(0),M2(0))")
#     
#     for (param in list(c(t0 = 0, t_jo = 0, t_1 = NA, t_2 = NA), 
#                        c(t0 = 1, t_jo = 0, t_1 = NA, t_2 = NA), 
#                        c(t0 = 1, t_jo = 1, t_1 = NA, t_2 = NA),
#                        c(t0 = 1, t_jo = NA, t_1 = 1, t_2 = 1),
#                        c(t0 = 1, t_jo = NA, t_1 = 0, t_2 = 1),
#                        c(t0 = 1, t_jo = NA, t_1 = 0, t_2 = 0)
#     )) {
#         t0 <- param["t0"]
#         if (!is.na(param["t_jo"])) {
#             t_jo <- param["t_jo"]
#             
#             # EIF calculation
#             ipwy <- (dat$A == t0) / g(t0, w)
#             hmjo <- pmaw(dat$M1, dat$M2, t_jo, w) / pm(dat$M1, dat$M2, dat$Z1, t0, w)
#             eify <- ipwy * hmjo / mean(ipwy * hmjo) * (dat$Y - my(dat$M1, dat$M2, dat$Z1, t0, w))
#             
#             ipwz <- ipwy
#             eifz <- ipwz / mean(ipwz) * (u(dat$Z, w, t0, t_jo) - intu(w, t0, t_jo))
#             
#             ipwm <- (dat$A == t_jo) / g(t_jo, w)
#             vbar <- intv(1, 1, w, t0) * pmaw(1, 1, t_jo, w) + intv(1, 0, w, t0) * pmaw(1, 0, t_jo, w) +
#                 intv(0, 1, w, t0) * pmaw(0, 1, t_jo, w) + intv(0, 0, w, t0) * pmaw(0, 0, t_jo, w)
#             eifm <- ipwm / mean(ipwm) * (intv(dat$M1, dat$M2, w, t0) - vbar)
#             
#             # eif <- rescale_y(eify + eifz + eifm + mubarZ[, 1], bounds$bounds)
#             eif <- eify + eifz + eifm + vbar
#             theta <- mean(eif)
#             
#             vname <- paste0("Y(", t0 ,",M(t_jo=", t_jo, "))")
#             vvbar[[vname]] <- mubarZ[, 1]
#             
#             thetas[vname] <- theta
#             eifs[[vname]] <- eif
#             
#         }
#         if (is.na(param[["t_jo"]])) {
#             t_1 <- param["t_1"]
#             t_2 <- param["t_2"]
#             
#             
#             
#             vname <- paste0("Y(", t0 ,",M1(", t_1, "),", "M2(", t_2, ")", ")")
#             # true_thetas[vname] <- weighted.mean(true_v, prob_w)
#             true_thetas[vname] <- mean(true_v)
#         }
#     }
#     
# }
# 
