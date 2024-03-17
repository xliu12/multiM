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

source("_simulation/fun_gendata.R")
devtools::load_all("MediationClustered")


GenData <- function(
        seedone = 123,
        num_clust = 100,
        clust_size = 30,
        quadratic.A = FALSE,
        quadratic.M = FALSE,
        quadratic.Y = FALSE,
        icca = 0.2,
        iccm = 0.2, 
        iccy = 0.2,
        x_z = 0,
        num_m = 2,
        num_x = 3,
        Yfamily = "gaussian"
) {

    N <- num_clust * clust_size

    data <- data.frame(
        id = 1:N,
        school = rep(1:num_clust, each = clust_size)
    )
    # Z (unobserved) cluster-level confounder -----------------
    z <- rep(rnorm(num_clust, sd = 1), each = clust_size)
    data$Z <- z
    
    # X individual-level confounders ------
    iccx <- 0.2
    gen_x <- list(iccx = iccx, x_z = x_z)
    
    xb <- mvtnorm::rmvnorm(n = num_clust,
                           mean = rep(0, num_x),
                           sigma = diag((1 - gen_x[["x_z"]]^2)*gen_x[["iccx"]], nrow = num_x))[rep(1:num_clust, each = clust_size), ]
    xe <- mvtnorm::rmvnorm(n = N,
                           mean = rep(0, num_x),
                           sigma = diag(1 - gen_x[["iccx"]], nrow = num_x))
    x <- gen_x[["x_z"]] * data$Z + xb + xe
    data$X <- x

    # Treatment ---------
    gen_a <- list(icca = icca, a_x = sqrt(0.15 * 1 / num_x), a_z = sqrt(0.4 / 1))
    ab <- rep(rnorm(num_clust, sd = sqrt(gen_a[["icca"]])), each = clust_size)
    if (quadratic.A == FALSE) {
        a_given <- ab +  gen_a[["a_x"]] * rowSums(data$X) + gen_a[["a_z"]] * data$Z
    }
    if (quadratic.A == TRUE) {
        Xquad <- (data$X ^ 2 - 1) / sqrt(2) # mean(rnorm(.)^2) = 1; var(rnorm(.)^2) = 2
        Xquad <- Xquad * data$Z
        gen_a[["a_x"]] <- sqrt(0.8 * 1 / ncol(Xquad))
        a_given <- ab +  gen_a[["a_x"]] * rowSums(Xquad) + gen_a[["a_z"]] * data$Z
    }
    data$A <- rbinom(N, 1, pa(1, a_given, gen_a))

    # Mediators ------------------
    iccm1 <- iccm2 <- iccm
    gen_m <- list(iccm1 = iccm1, iccm2 = iccm2,
                  m1_on_a = 0.2, m1_on_x = sqrt(0.15 / num_x), m1_on_z = sqrt(0.4),
                  m2_on_a = 0.2, m2_on_m1 = 0.2, 
                  m2_on_am1 = 0.2,
                  m2_on_x = sqrt(0.15 / num_x),
                  m2_on_z = sqrt(0.4)
    )

    m1b <- rep(rnorm(num_clust, sd = sqrt(gen_m[["iccm1"]])), each = clust_size)
    m2b <- rep(rnorm(num_clust, sd = sqrt(gen_m[["iccm2"]])), each = clust_size)

    if (quadratic.M == TRUE) {
        Xquad <- (data$X ^ 2 - 1) / sqrt(2)
        Xquad <- Xquad * data$Z
        gen_m[["m1_on_x"]] <- sqrt(0.8 * 1 / ncol(Xquad))
        gen_m[["m2_on_x"]] <- sqrt(0.8 * 1 / ncol(Xquad))
        
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
        data$M1 <- rbinom(N, 1, prob = pm1(1, data$A, m1_given, gen_m) )

        m2_given <- m2b + gen_m[["m2_on_x"]] * rowSums(data$X) +
            gen_m[["m2_on_z"]] * data$Z
        data$M2 <- rbinom(N, 1, prob = pm2(1, data$M1, data$A, m2_given, gen_m) )
    }

    # Outcome -------------------

    y_on_amint <- 0
    gen_y <- list(iccy = iccy, yintercept = 1,
                  y_on_a = 0.2, y_on_m1 = 0.2, y_on_m2 = 0.2,
                  y_on_am1 = y_on_amint, y_on_am2 = y_on_amint, 
                  y_on_m1m2 = 0.2, y_on_am1m2 = 0.2,
                  y_on_x = sqrt(0.15 / num_x), y_on_z = sqrt(0.4)
                  )
    yb <- rnorm(num_clust, sd = sqrt(gen_y[["iccy"]]))[rep(1:num_clust, each = clust_size)]

    if (quadratic.Y == TRUE) {
        Xquad <- (data$X ^ 2 - 1) / sqrt(2)
        Xquad <- Xquad * data$Z
        gen_y[["y_on_x"]] <- sqrt(0.8 * 1 / ncol(Xquad))

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






# simulation conditions ----------------------------

condition_all <- data.frame(expand.grid(
    num_clust = c(100, 20),
    clust_size = c(10, 20, 50),
    # generating data
    quadratic.A = c(F, T),
    quadratic.M = c(F, T),
    quadratic.Y = c(F, T),
    icc = c(0.2, 0.4),
    x_z = c(0, 0.3),
    Yfamily = c("gaussian") 
))

# specific conditions for a job
condition <- condition_all %>%
    filter(
        quadratic.A + quadratic.M + quadratic.Y == 0,
        x_z == 0,
        icc == 0.2, num_clust == 100, clust_size %in% c(20)
        ) 

methds_all <- data.frame(expand.grid(
    cluster_a = c("FE", "RE", "noncluster"), # fixed-effects, random-effects, single-level models for the treatment
    cluster_m =  c("FE", "RE", "noncluster"), # fixed-effects, random-effects, single-level models for the mediators
    cluster_y =  c("FE", "RE", "noncluster"), # fixed-effects, random-effects, single-level models for the outcome
    Fit = c("glm", "mlr") # model estimation: parametric model with linear terms, nonparametric with superlearner package 
    ))
# methods considered in a job
methds <- methds_all %>% 
    filter(Fit == "glm") %>%
    mutate(
        cluster_opt_a = paste(cluster_a, Fit, sep="."),
        cluster_opt_m = paste(cluster_m, Fit, sep="."),
        cluster_opt_y = paste(cluster_y, Fit, sep=".")
    ) 

set.seed(12)
datseeds <- sample(1:1e6, 1000)

iseed <- 1
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
        icca = condition$icc[cond],
        iccm = condition$icc[cond],
        iccy = condition$icc[cond],
        x_z = condition$x_z[cond],
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

    # true values to be averaged over across 1000 replications of the large sample size
    gen_largeN <- GenData(
        seedone = datseeds[iseed],
        num_clust = 500,
        clust_size = 500,
        quadratic.A = condition$quadratic.A[cond],
        quadratic.M = condition$quadratic.M[cond],
        quadratic.Y = condition$quadratic.Y[cond],
        icca = condition$icc[cond],
        iccm = condition$icc[cond],
        iccy = condition$icc[cond],
        x_z = condition$x_z[cond],
        Yfamily = Yfamily,
        num_m = 2,
        num_x = 3
    )
    true_values <- unlist(gen_largeN$truevals)

    # estimators 
    one.jj <- function(jj = 1) {
        Fit <- as.character(methds$Fit[jj])
        if (Fit == "glm") {
            learners_a <- learners_m <- learners_y <- c("SL.glm")
        }
        if (Fit == "mlr") {
            learners_a <- learners_m <- learners_y <- c("SL.ranger", "SL.earth")
        }
        
        Wnames <- NULL
        cluster_opt_a <- methds$cluster_opt_a[jj]
        cluster_opt_m <- methds$cluster_opt_m[jj]
        cluster_opt_y <- methds$cluster_opt_y[jj]
        
        out <- clustered(data,
                         Sname,
                         Wnames, 
                         Xnames,
                         Aname,
                         Mnames, 
                         Yname, Yfamily,
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
            est_multi = out$thetas_effs, # multiply-robust estimator
            est_reg = out$regs_effs, # regression-based estimator
            est_rmpw = out$rmpw_effs # weighting-based estimator
            )
        estimates1 <- estimates[str_detect(rownames(estimates), "Y"), ]
        
        res <- data.frame(methds[jj, ],
                          estimand = names(true_values),
                          true_val = true_values,
                          estimates1, row.names = NULL)
        
        return(res)
    }
    
    oneboot <- lapply(1:nrow(methds), one.jj)
    one_boot <- do.call(rbind, oneboot)
    
    res <- data.frame(condition[cond, ],
                      one_boot, 
                      row.names = NULL)

    return(res)
}


# running ----

set.seed(12)
jobconds <- c(1:nrow(condition))

reslist <- list()

# try a few replications
# jobseeds <- 1:10
jobseeds <- 1:1000

nj <- paste(unique(condition$clust_size), collapse = "-")
icc <- paste(unique(condition$icc)*10, collapse = "-")
J <- paste(unique(condition$num_clust), collapse = "-")
learner <- unique(methds$Fit)
xz <- paste(unique(condition$x_z)*10, collapse = "-")
rname <- glue("_simulation/xz{xz}_Fit_{learner}_icc{icc}_J{J}_n{nj}_rep{jobseeds[1]}_{tail(jobseeds, 1)}.RData")
rname <- as.character(rname)

for(k in jobconds) {
    res <- mclapply(jobseeds, 
                    OneData, cond = k, 
                    mc.preschedule = TRUE, mc.cores = 100) 
    reslist[[k]] <- res
    save(reslist, file = as.character(rname))
}


# resulting ----
dirpath_s <- rname
read.one <- function(dirpath) {
    load(file.path(dirpath))
    res_bind <- lapply(
        reslist[sapply(reslist, length) > 0], 
        function(l){
            do.call(rbind, l[sapply(l, length) > 1])
        })
    
    resdf <- do.call(rbind, res_bind)
    resdf1 <- resdf
    return(resdf1)
}
resdf <- lapply(dirpath_s, read.one)
resdf1 <- do.call(rbind, resdf)

res_df <- resdf1 %>%
    group_by(quadratic.A, quadratic.M, quadratic.Y, icc, Yfamily, x_z, estimand) %>% 
    mutate(
        trueval = mean(true_val) 
    ) %>% 
    group_by(num_clust, clust_size, quadratic.A, quadratic.M, quadratic.Y, icc, Yfamily, x_z, estimand, cluster_a, cluster_m, cluster_y, Fit) %>% 
    mutate(
        rbias = across(contains("est_"), ~. / trueval - 1),
        bias = across(contains("est_"), ~. - trueval),
        MSE = across(contains("est_"), ~(. - trueval)^2),
        nreps = n()
    ) 

res_df <- do.call(data.frame, res_df)


res_summ <- res_df %>%
    group_by(num_clust, clust_size, quadratic.A, quadratic.M, quadratic.Y, icc, Yfamily, x_z, estimand, cluster_a, cluster_m, cluster_y, Fit) %>% 
    summarise(across(c(contains("bias"), contains("MSE"), nreps), ~mean(.)))  

res_summ <- res_summ %>% mutate(
    across(contains("MSE"), .fns = list(sqrt = ~sqrt(.)), 
    .names ="R{.col}")
)
# _reg aggregate over a_model (a_model is not involved in the regression-based estimator)
res_summ1 <- res_summ %>% 
    group_by(num_clust, clust_size, quadratic.A, quadratic.M, quadratic.Y, icc, Yfamily, x_z, estimand, cluster_m, cluster_y, Fit # cluster_a 
             ) %>% 
    mutate(across(ends_with("_reg"), ~mean(.)))
# _rmpw aggregate over y_model (y_model is not involved in the weighting-based estimator)
res_summ1 <- res_summ1 %>% 
    ungroup() %>% 
    group_by(num_clust, clust_size, quadratic.A, quadratic.M, quadratic.Y, icc, Yfamily, x_z, estimand, cluster_m,  cluster_a, Fit #cluster_y
             ) %>% 
    mutate(across(ends_with("_rmpw"), ~mean(.)))

res_long <- res_summ1 %>% 
    pivot_longer(
        cols = c(ends_with("_multi"), ends_with("_reg"), ends_with("_rmpw")), 
        names_to = c("performance", "estimator"),
        names_pattern = "(.*).est_(.*)"
    )
res_long1 <- res_long %>% mutate(
    performance = factor(performance),
    J = factor(num_clust), nj = factor(clust_size),
    A_model = case_match(cluster_a,
                           "noncluster" ~ "SL", .default = cluster_a),
    M_model = case_match(cluster_m,
                           "noncluster" ~ "SL", .default = cluster_m),
    Y_model = case_match(cluster_y,
                           "noncluster" ~ "SL", .default = cluster_y),
    A_mod = ifelse(cluster_a %in% c("FE","RE"), "ML", "SL"),
    M_mod = ifelse(cluster_m %in% c("FE","RE"), "ML", "SL"),
    Y_mod = ifelse(cluster_y %in% c("FE","RE"), "ML", "SL"),
    abs.value = abs(value),
    estimand = case_match(
        estimand,
        "Y(0,gmjo(0))" ~ "theta_joint(0,0)",
        "Y(1,gmjo(0))" ~ "theta_joint(1,0)",
        "Y(1,gmjo(1))" ~ "theta_joint(1,1)",
        "Y(1,gm1(0),gm2(0))" ~ "theta_each(1,0,0)",
        "Y(1,gm1(1),gm2(0))" ~ "theta_each(1,1,0)",
        "Y(1,gm1(1),gm2(1))" ~ "theta_each(1,1,1)"
        ),
    estimator = case_match(
        estimator,
        "multi" ~ "multiply_robust", 
        "reg" ~ "regression_based",
        "rmpw" ~ "weighting_based"
    ),
    abs.performance = case_when(
        as.character(performance)=="bias" ~ "|Bias|",
        as.character(performance)=="rbias" ~ "|Rel.Bias|",
        as.character(performance)=="RMSE" ~ "RMSE",
        as.character(performance)=="MSE" ~ "MSE"
    )
) %>% 
    mutate(abs.performance = factor(abs.performance, levels = c("|Rel.Bias|", "RMSE", "|Bias|", "MSE")))

