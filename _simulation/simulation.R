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
library(nnet)
library(gam)
library(lightgbm)
library(origami)




GenData <- function(
        seedone = 123,
        J = 100,
        nj = 30,
        quadratic.A = FALSE,
        quadratic.M = FALSE,
        quadratic.Y = FALSE,
        if.null = FALSE,
        icca = 0.2,
        iccm = 0.2, 
        iccy = 0.2,
        x_z = 0,
        num_m = 2,
        num_x = 3,
        Yfamily = "gaussian",
        indep_M = FALSE,
        int.XZ = TRUE
) {
    
    set.seed(seed = seedone)
    
    if (nj==5) { 
        njrange <- c(5, 20)
    }
    
    
    if (nj == 50) {
        njrange <- c(50, 100)
    }
    
    nj_sizes <- runif(J, njrange[1], njrange[2]) %>% round()
    
    N <- sum(nj_sizes)
    
    data <- data.frame(
        id = 1:N,
        school = unlist(map(1:J, ~rep(.x, each = nj_sizes[.x])))
    ) %>% 
        group_by(school) %>% 
        mutate(W_nj = (n()-njrange[1])/(njrange[2]-njrange[1]) )
    
    
    
    # Z (unobserved) cluster-level confounders -----------------
    data$Z <- unlist(map(1:J, ~rep(rnorm(1), each = nj_sizes[.x])))
    
    # X individual-level confounders ------
    iccx <- 0.2
    gen_x <- list(iccx = iccx, x_z = x_z)
    
    xb <- mvtnorm::rmvnorm(n = J,
                           mean = rep(0, num_x),
                           sigma = diag((1 - gen_x[["x_z"]]^2)*gen_x[["iccx"]], nrow = num_x))[data$school, ] #[rep(1:J, each = nj), ]
    xe <- mvtnorm::rmvnorm(n = N,
                           mean = rep(0, num_x),
                           sigma = diag(1 - gen_x[["iccx"]], nrow = num_x))
    x <- gen_x[["x_z"]] * data$Z + xb + xe
    data$X <- x
    
    
    # Treatment ---------
    # icca <- 0.2
    gen_a <- list(icca = icca, a_x = sqrt(0.15 * 1 / num_x), a_z = sqrt(0.4 / 1))
    ab <- unlist(map(1:J, ~rep(rnorm(1, mean = 0, sd = sqrt(gen_a[["icca"]])), each = nj_sizes[.x])))
    if (quadratic.A == FALSE) {
        Xlinear <- data$X
        
        a_given <- ab +  gen_a[["a_x"]] * rowSums(Xlinear) + gen_a[["a_z"]] * data$Z
    }
    if (quadratic.A == TRUE) {
        Xquad <- (data$X ^ 2 - 1) / sqrt(4) 
        a_given <- ab +  gen_a[["a_x"]] * rowSums(Xquad) + gen_a[["a_z"]] * data$Z
    }
    
    data$A <- rbinom(N, 1, pa(1, a_given, gen_a))
    
    # Mediators ------------------
    # iccm <- 0.2
    iccm1 <- iccm2 <- iccm
    gen_m <- list(iccm1 = iccm1, iccm2 = iccm2,
                  m1_on_a = 0.2, m1_on_x = sqrt(0.15 / num_x), m1_on_z = sqrt(0.4),
                  m1_on_az = 0.2, m1_on_anj = 0.2,
                  # M2
                  m2_on_a = 0.2, m2_on_m1 = 0.2, m2_on_am1 = 0.2,
                  m2_on_x = sqrt(0.15 / num_x), m2_on_z = sqrt(0.4),
                  m2_on_az = 0.2, m2_on_anj = 0.2
    )
    # conditional independent mediators
    if (indep_M) {
        gen_m[c("m2_on_m1","m2_on_am1")] <- 0
    }
    
    if (int.XZ==FALSE) {
        # gen_m[c("m1_on_az","m2_on_az")] <- 0
        gen_m[c("m1_on_anj","m2_on_anj")] <- 0
    }
    m1b <- unlist(map(1:J, ~rep(rnorm(1, mean = 0, sd = sqrt(gen_m[["iccm1"]])), each = nj_sizes[.x])))
    m2b <- unlist(map(1:J, ~rep(rnorm(1, mean = 0, sd = sqrt(gen_m[["iccm2"]])), each = nj_sizes[.x])))
    
    if (quadratic.M == TRUE) {
        Xquad <- (data$X ^ 2 - 1) / sqrt(4)
        
        m1_given <- m1b + gen_m[["m1_on_x"]] * rowSums(Xquad) +
            gen_m[["m1_on_z"]] * data$Z
        data$M1 <- rbinom(N, 1, prob = pm1(1, data$A, data$Z, data$W_nj, m1_given, gen_m) )
        
        m2_given <- m2b + gen_m[["m2_on_x"]] * rowSums(Xquad) +
            gen_m[["m2_on_z"]] * data$Z
        data$M2 <- rbinom(N, 1, prob = pm2(1, data$M1, data$A, data$Z, data$W_nj, m2_given, gen_m) )
        
    }
    if (quadratic.M == FALSE) {
        Xlinear <- data$X
        
        m1_given <- m1b + gen_m[["m1_on_x"]] * rowSums(Xlinear) +
            gen_m[["m1_on_z"]] * data$Z
        data$M1 <- rbinom(N, 1, prob = pm1(1, data$A, data$Z, data$W_nj, m1_given, gen_m) )
        
        m2_given <- m2b + gen_m[["m2_on_x"]] * rowSums(Xlinear) +
            gen_m[["m2_on_z"]] * data$Z
        data$M2 <- rbinom(N, 1, prob = pm2(1, data$M1, data$A, data$Z, data$W_nj, m2_given, gen_m) )
    }
    
    # Outcome -------------------
    
    y_on_amint <- 0
    # iccy <- 0.2
    gen_y <- list(iccy = iccy, yintercept = 1,
                  y_on_a = 0.2, y_on_m1 = 1, y_on_m2 = 1,
                  y_on_am1 = 0, y_on_am2 = 0, 
                  y_on_m1m2 = 0.2, y_on_am1m2 = 0.2,
                  y_on_az = 0.2, y_on_m1z = 0.2, y_on_m2z = 0.2,
                  y_on_anj = 0.2,  
                  y_on_x = sqrt(0.15 / num_x), y_on_z = sqrt(0.4)
    )
    yb <- rnorm(J, sd = sqrt(gen_y[["iccy"]]))[data$school] #[rep(1:J, each = nj)]
    
    if (if.null==TRUE) {
        gen_y <- list(iccy = iccy, yintercept = 1,
                      # y_on_a = 0.2, y_on_m1 = 0.2, y_on_m2 = 0.2,
                      y_on_a = 0, y_on_m1 = 0, y_on_m2 = 0,
                      y_on_am1 = 0, y_on_am2 = 0, 
                      y_on_m1m2 = 0, y_on_am1m2 = 0,
                      # y_on_m1m2 = 0, y_on_am1m2 = 0,
                      y_on_az = 0, y_on_m1z = 0, y_on_m2z = 0,
                      y_on_anj = 0,  
                      y_on_x = sqrt(0.15 / num_x), y_on_z = sqrt(0.4)
        )
    }
    
    if (int.XZ==FALSE) {
        gen_y[c("y_on_anj")] <- 0
    }
    
    if (quadratic.Y == TRUE) {
        Xquad <- (data$X ^ 2 - 1) / sqrt(4)
        
        y_given <- yb + gen_y[["yintercept"]] + gen_y[["y_on_x"]] * rowSums(Xquad) +
            gen_y[["y_on_z"]] * data$Z
    }
    if (quadratic.Y == FALSE) {
        Xlinear <- data$X
        
        y_given <- yb + gen_y[["yintercept"]] + gen_y[["y_on_x"]] * rowSums(Xlinear) +
            gen_y[["y_on_z"]] * data$Z
    }
    
    if (Yfamily == "gaussian") {
        condmy <- my(m2 = data$M2, m1 = data$M1, data$A, data$Z, data$W_nj, given = y_given, gen_y, binary = FALSE)
        data$Y <- condmy + rnorm(N, sd = sqrt(1 - gen_y[["iccy"]]))
    }
    if (Yfamily == "binomial") {
        condmy <- my(m2 = data$M2, m1 = data$M1, data$A, data$Z, data$W_nj, given = y_given, gen_y, binary = TRUE)
        data$Y <- rbinom(N, 1, condmy)
    }
    
    datobs <- do.call(data.frame, data)
    
    
    # nj_sizes <- as.numeric(table(data$school))
    
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
        truevals <- truevals_individual <- truevals_cluster <- list()
        j <- 1
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
                    p_marginal <- pm1(m1, a1, data$Z, data$W_nj, m1_given, gen_m) *
                        pm2a(m2, a2, data$Z, data$W_nj, m2_given, m1_given, gen_m)
                    
                    p_marginal * my(m2 = vals$m2[i], m1 = vals$m1[i], a = a0, data$Z, data$W_nj,
                                    y_given, gen_y, binary = (Yfamily=="binomial"))
                    
                }, simplify = TRUE)
                int_m1m2 <- rowSums(intm1m2)
                truevals_individual[[glue("Y({a0},gm1({a1}),gm2({a2}))")]] <- mean(int_m1m2)
                truevals_cluster[[glue("Y({a0},gm1({a1}),gm2({a2}))")]] <- mean(aggregate(data.frame(int_m1m2), by=list(data$school), mean)[,2])
            }
            
            if (!is.na(ajo)) {
                # integrate over the joint mediator distribution, (M1, M2) under ajo
                intm1m2 <- sapply(1:nrow(vals), function(i=1) {
                    m1 <- vals$m1[i]
                    m2 <- vals$m2[i]
                    p_joint <- pm1(m1, ajo, data$Z, data$W_nj, m1_given, gen_m) *
                        pm2(m2, m1, ajo, data$Z, data$W_nj, m2_given, gen_m)
                    
                    p_joint * my(m2 = vals$m2[i], m1 = vals$m1[i], a = a0, data$Z, data$W_nj,
                                 y_given, gen_y, binary = (Yfamily=="binomial"))
                    
                }, simplify = TRUE)
                int_m1m2 <- rowSums(intm1m2)
                truevals_individual[[glue("Y({a0},gmjo({ajo}))")]] <- mean(int_m1m2)
                truevals_cluster[[glue("Y({a0},gmjo({ajo}))")]] <- mean(aggregate(data.frame(int_m1m2), by=list(data$school), mean)[,2])
            }
        }
        
        
        list(truevals_individual=truevals_individual, truevals_cluster=truevals_cluster)
    }
    
    truevals <- trueVals()
    
    out <- mget(ls(envir = environment()))
    
    return(out)
    
}






# Simulation----------------------------

condition_all <- data.frame(expand.grid(
    J = c(100, 70,40,  20, 10),
    nj = c(5),
    # generating data
    quadratic.A = c(F, T),
    quadratic.M = c(F, T),
    quadratic.Y = c(F, T),
    if.null = c(F,T)
))


condition <- condition_all %>%
    filter((quadratic.A + quadratic.M + quadratic.Y == 0) | (quadratic.A + quadratic.M + quadratic.Y == 3)) %>%
    filter( if.null == F) %>% 
    filter(J %in% c(100, 70,40)
           # J %in% c(40,  20, 10) 
    ) %>% 
    mutate(nj = 5 ) # nj=50


methds_all <- data.frame(expand.grid(
    Morder = c("21", "12"),
    Fit = c("mlr","glm"), 
    cluster_opt = c("cwc.FE", "cwc") 
)) 
methds <- methds_all 


datseeds <- c(sample(1:1e6, 1000))

iseed <-1
cond <- 1

source("_simulation/fun_gen.R")

OneData <- function(iseed = 1, cond = 1){
    seedone <- datseeds[iseed]
    gen_data <- GenData(
        seedone = seedone, # datseeds[iseed],
        J = condition$J[cond],
        nj = condition$nj[cond],
        quadratic.A = condition$quadratic.A[cond],
        quadratic.M = condition$quadratic.M[cond],
        quadratic.Y = condition$quadratic.Y[cond],
        if.null = condition$if.null[cond]
    )
    
    data <- gen_data$datobs
    
    
    Sname <- "school"
    Xnames <- colnames(data)[grep("^X", colnames(data))]
    Wnames <- colnames(data)[grep("^W", colnames(data))]
    
    Aname <- "A"
    Mnames <- colnames(data)[grep("^M", colnames(data))]
    Yname <- "Y"
    
    
    
    # true_values <- unlist(gen_largeJ$truevals)
    # trueEffs_individual <- effect(gen_largeJ$truevals$truevals_individual) %>% unlist()
    # trueEffs_cluster <- effect(gen_largeJ$truevals$truevals_cluster) %>% unlist()
    
    # estimators -----
    jj <- 1
    one.jj <- function(jj = 1) {
        Fit <- as.character(methds$Fit[jj])
        if (Fit == "glm") {
            learners_a <- learners_m <- learners_y <- c("SL.glm")
            num_folds <- 1
        }
        if (Fit == "mlr") {
            learners_a <- learners_m <- learners_y <- c("SL.nnet", "SL.gam")
            
            num_folds <- 5 
        }
        
        
        out <- ClusterMed(data = data,
                          Sname = Sname,
                          Wnames = Wnames, 
                          Xnames = Xnames,
                          Aname = Aname,
                          Mnames = Mnames, 
                          Yname = Yname, 
                          cluster_opt = methd$cluster_opt[jj],
                          Morder = methd$Morder[jj],
                          num_folds = num_folds,
                          learners_a = learners_a,  
                          learners_m = learners_m,
                          learners_y = learners_y )
        
        
        estimates <- out
        
        res <- data.frame(methds[jj, ],
                          estimates, row.names = NULL)
        
        
        
        return(res)
    }
    
    oneboot <- lapply(1:nrow(methds), one.jj)
    one_boot <- do.call(rbind, oneboot)
    
    res <- data.frame(condition[cond, ],
                      one_boot, 
                      row.names = NULL)
    
    
    return(res)
}

