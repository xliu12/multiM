mu <- function(data, varnames, family, folds, learners, ...) {
    b <- matrix(nrow = nrow(data), ncol = 2)
    colnames(b) <- c("mu(0,Z,M1,M2,C)", "mu(1,Z,M1,M2,C)")
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid_1 <- valid_0 <- origami::validation(data, folds[[v]])
        valid_1[[varnames$A]] <- 1
        valid_0[[varnames$A]] <- 0
        
        preds <- crossfit(train, list(valid_0, valid_1), varnames$Y,
                          c(varnames$W, varnames$A, varnames$Z, varnames$M1, varnames$M2),
                          family, learners = learners)
        b[folds[[v]]$validation_set, 1] <- preds[[1]]
        b[folds[[v]]$validation_set, 2] <- preds[[2]]
    }
    b
}



# pseudo outcome
# integrating over marginal mediator distributions M1 * M2

mu_M1M2 <- function(data, varnames, bb, hm1m2, t0, folds, learners, ...) {
    u <- matrix(nrow = nrow(data), ncol = 1)
    colnames(u) <- "mu_M1M2(Z,C)"
    
    data[["mu(t0,Z,M1,M2,C)h_M1M2"]] <- bb[, gl("mu({t0},Z,M1,M2,C)")]*hm1m2
    
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[varnames$A]] <- t0
        
        u[folds[[v]]$validation_set, "mu_M1M2(Z,C)"] <-
            crossfit(train, list(valid), "mu(t0,Z,M1,M2,C)h_M1M2",
                     c(varnames$Z, varnames$A, varnames$W),
                     "gaussian", learners = learners)[[1]]
    }
    u
}

mubar_M1M2 <- function(data, varnames, muM1M2, t0, folds, learners, ...) {
    ubar <- matrix(nrow = nrow(data), ncol = 1)
    colnames(ubar) <- "mubar_M1M2(C)"
    
    data[["mu_M1M2(Z,C)"]] <- muM1M2[, "mu_M1M2(Z,C)"]
    
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[varnames$A]] <- t0
        
        ubar[folds[[v]]$validation_set, "mubar_M1M2(C)"] <-
            crossfit(train, list(valid), "mu_M1M2(Z,C)",
                     c(varnames$A, varnames$W),
                     "gaussian", learners = learners)[[1]]
    }
    ubar
}


# integrating over p(Z|t0,C) * p(M1|t_1, C) 
mu_ZM1 <- function(data, varnames, bb, hm1z, t0, folds, learners, ...) {
    u <- matrix(nrow = nrow(data), ncol = 1)
    colnames(u) <- "mu_ZM1(M2,C)"
    
    data[["mu(t0,Z,M1,M2,C)h_M1Z"]] <- bb[, gl("mu({t0},Z,M1,M2,C)")]*hm1z
    
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[varnames$A]] <- t0
        
        u[folds[[v]]$validation_set, "mu_ZM1(M2,C)"] <-
            crossfit(train, list(valid), "mu(t0,Z,M1,M2,C)h_M1Z",
                     c(varnames$M2, varnames$A, varnames$W),
                     "gaussian", learners = learners)[[1]]
    }
    u
}

mubar_ZM1 <- function(data, varnames, muZM1, t_2, folds, learners, ...) {
    ubar <- matrix(nrow = nrow(data), ncol = 1)
    colnames(ubar) <- "mubar_ZM1(C)"
    
    data[["mu_ZM1(M2,C)"]] <- muZM1[, "mu_ZM1(M2,C)"]
    
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[varnames$A]] <- t_2
        
        ubar[folds[[v]]$validation_set, "mubar_ZM1(C)"] <-
            crossfit(train, list(valid), "mu_ZM1(M2,C)",
                     c(varnames$A, varnames$W),
                     "gaussian", learners = learners)[[1]]
    }
    ubar
}


# integrating over p(Z|t0,C) * p(M2|t_2, C) 
mu_ZM2 <- function(data, varnames, bb, hm2z, t0, folds, learners, ...) {
    u <- matrix(nrow = nrow(data), ncol = 1)
    colnames(u) <- "mu_ZM2(M1,C)"
    
    data[["mu(t0,Z,M1,M2,C)h_M2Z"]] <- bb[, gl("mu({t0},Z,M1,M2,C)")]*hm2z
    
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[varnames$A]] <- t0
        
        u[folds[[v]]$validation_set, "mu_ZM2(M1,C)"] <-
            crossfit(train, list(valid), "mu(t0,Z,M1,M2,C)h_M2Z",
                     c(varnames$M1, varnames$A, varnames$W),
                     "gaussian", learners = learners)[[1]]
    }
    u
}

mubar_ZM2 <- function(data, varnames, muZM2, t_1, folds, learners, ...) {
    ubar <- matrix(nrow = nrow(data), ncol = 1)
    colnames(ubar) <- "mubar_ZM2(C)"
    
    data[["mu_ZM2(M1,C)"]] <- muZM2[, "mu_ZM2(M1,C)"]
    
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[varnames$A]] <- t_1
        
        ubar[folds[[v]]$validation_set, "mubar_ZM2(C)"] <-
            crossfit(train, list(valid), "mu_ZM2(M1,C)",
                     c(varnames$A, varnames$W),
                     "gaussian", learners = learners)[[1]]
    }
    ubar
}


# using observed Z distribution
mu_Mjo <- function(data, varnames, bb, hmjo, t0, folds, learners, ...) {
    u <- matrix(nrow = nrow(data), ncol = 1)
    colnames(u) <- "mu_Mjo(Z,C)"
    
    data[["mu(t0,Z,M1,M2,C)h_Mjo"]] <- bb[, gl("mu({t0},Z,M1,M2,C)")]*hmjo
    
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[varnames$A]] <- t0

        u[folds[[v]]$validation_set, "mu_Mjo(Z,C)"] <-
            crossfit(train, list(valid), "mu(t0,Z,M1,M2,C)h_Mjo",
                     c(varnames$Z, varnames$A, varnames$W),
                     "gaussian", learners = learners)[[1]]
    }
    u
}

mubar_Mjo <- function(data, varnames, muMjo, t0, folds, learners, ...) {
    ubar <- matrix(nrow = nrow(data), ncol = 1)
    colnames(ubar) <- "mubar_Mjo(C)"
    
    data[["mu_Mjo(Z,C)"]] <- muMjo[, "mu_Mjo(Z,C)"]
    
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[varnames$A]] <- t0

        ubar[folds[[v]]$validation_set, "mubar_Mjo(C)"] <-
            crossfit(train, list(valid), "mu_Mjo(Z,C)",
                     c(varnames$A, varnames$W),
                     "gaussian", learners = learners)[[1]]
    }
    ubar
}


# using observed joint mediator distribution
mu_Z <- function(data, varnames, bb, hz, t0, folds, learners, ...) {
    vv <- matrix(nrow = nrow(data), ncol = 1)
    colnames(vv) <- "mu_Z(Mjo,C)"
    
    data[["mu(t0,Z,M1,M2,C)hz"]] <- bb[, gl("mu({t0},Z,M1,M2,C)")]*hz[, gl("h_z({t0})")]
    
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[varnames$A]] <- t0

        vv[folds[[v]]$validation_set, "mu_Z(Mjo,C)"] <-
            crossfit(train, list(valid), "mu(t0,Z,M1,M2,C)hz",
                     c(varnames$M1, varnames$M2, varnames$A, varnames$W), "gaussian",
                     learners = learners)[[1]]
    }
    vv
}

mubar_Z <- function(data, varnames, muZ, t_jo, folds, learners, ...) {
    vbar <- matrix(nrow = nrow(data), ncol = 1)
    colnames(vbar) <- "mubar_Z(C)"
    
    data[["mu_Z(Mjo,C)"]] <- muZ[, "mu_Z(Mjo,C)"]
    
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[varnames$A]] <- t_jo

        vbar[folds[[v]]$validation_set, "mubar_Z(C)"] <-
            crossfit(train, list(valid), "mu_Z(Mjo,C)",
                     c(varnames$A, varnames$W), "gaussian",
                     learners = learners)[[1]]
    }
    vbar
}
