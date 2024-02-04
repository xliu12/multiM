h_zm1 <- function(data, varnames, folds, learners, ...) {
    h_zm1 <- matrix(nrow = nrow(data), ncol = 2)
    colnames(h_zm1) <- c("h_zm1(0)", "h_zm1(1)")
    for (v in seq_along(folds)) {
        # train <- stack_data(origami::training(data, folds[[v]]), npsem)
        train <- stack_data(origami::training(data, folds[[v]]), varnames$M2)
        valid_1 <- valid_0 <- origami::validation(data, folds[[v]])
        valid_1[[varnames$A]] <- 1
        valid_0[[varnames$A]] <- 0

        p_zm1m2w <- crossfit(train, list(valid_0, valid_1), "tmp_tmce_delta",
                          c(varnames$Z, varnames$M1, varnames$M2, varnames$W, varnames$A),
                          id = "tmp_tmce_id",
                          "binomial", learners = learners, bound = TRUE)

        p_m2w <- crossfit(train, list(valid_0, valid_1), "tmp_tmce_delta",
                         c(varnames$M2, varnames$W, varnames$A),
                         id = "tmp_tmce_id",
                         "binomial", learners = learners, bound = TRUE)

        h_zm1[folds[[v]]$validation_set, "h_zm1(0)"] <-
            (p_zm1m2w[[1]] / (1 - p_zm1m2w[[1]])) * ((1 - p_m2w[[1]]) / p_m2w[[1]])

        h_zm1[folds[[v]]$validation_set, "h_zm1(1)"] <-
            (p_zm1m2w[[2]] / (1 - p_zm1m2w[[2]])) * ((1 - p_m2w[[2]]) / p_m2w[[2]])
    }
    h_zm1
}

h_zm2 <- function(data, varnames, folds, learners, ...) {
    h_zm2 <- matrix(nrow = nrow(data), ncol = 2)
    colnames(h_zm2) <- c("h_zm2(0)", "h_zm2(1)")
    for (v in seq_along(folds)) {
        # train <- stack_data(origami::training(data, folds[[v]]), npsem)
        train <- stack_data(origami::training(data, folds[[v]]), varnames$M1)
        valid_1 <- valid_0 <- origami::validation(data, folds[[v]])
        valid_1[[varnames$A]] <- 1
        valid_0[[varnames$A]] <- 0
        
        p_zm1m2w <- crossfit(train, list(valid_0, valid_1), "tmp_tmce_delta",
                             c(varnames$Z, varnames$M1, varnames$M2, varnames$W, varnames$A),
                             id = "tmp_tmce_id",
                             "binomial", learners = learners, bound = TRUE)
        
        p_m1w <- crossfit(train, list(valid_0, valid_1), "tmp_tmce_delta",
                          c(varnames$M1, varnames$W, varnames$A),
                          id = "tmp_tmce_id",
                          "binomial", learners = learners, bound = TRUE)
        
        h_zm2[folds[[v]]$validation_set, "h_zm2(0)"] <-
            (p_zm1m2w[[1]] / (1 - p_zm1m2w[[1]])) * ((1 - p_m1w[[1]]) / p_m1w[[1]])
        
        h_zm2[folds[[v]]$validation_set, "h_zm2(1)"] <-
            (p_zm1m2w[[2]] / (1 - p_zm1m2w[[2]])) * ((1 - p_m1w[[2]]) / p_m1w[[2]])
    }
    h_zm2
}

h_z1 <- function(data, varnames, folds, learners, ...) {
    h_z1 <- matrix(nrow = nrow(data), ncol = 2)
    colnames(h_z1) <- c("h_z1(0)", "h_z1(1)")
    for (v in seq_along(folds)) {
        # train <- stack_data(origami::training(data, folds[[v]]), npsem)
        train <- stack_data(origami::training(data, folds[[v]]), varnames$M1)
        valid_1 <- valid_0 <- origami::validation(data, folds[[v]])
        valid_1[[varnames$A]] <- 1
        valid_0[[varnames$A]] <- 0
        
        p_zm1w <- crossfit(train, list(valid_0, valid_1), "tmp_tmce_delta",
                             c(varnames$Z, varnames$M1, varnames$W, varnames$A),
                             id = "tmp_tmce_id",
                             "binomial", learners = learners, bound = TRUE)
        
        p_m1w <- crossfit(train, list(valid_0, valid_1), "tmp_tmce_delta",
                          c(varnames$M1, varnames$W, varnames$A),
                          id = "tmp_tmce_id",
                          "binomial", learners = learners, bound = TRUE)
        
        h_z1[folds[[v]]$validation_set, "h_z1(0)"] <-
            (p_zm1w[[1]] / (1 - p_zm1w[[1]])) * ((1 - p_m1w[[1]]) / p_m1w[[1]])
        
        h_z1[folds[[v]]$validation_set, "h_z1(1)"] <-
            (p_zm1w[[2]] / (1 - p_zm1w[[2]])) * ((1 - p_m1w[[2]]) / p_m1w[[2]])
    }
    h_z1
}


h_z2 <- function(data, varnames, folds, learners, ...) {
    h_z2 <- matrix(nrow = nrow(data), ncol = 2)
    colnames(h_z2) <- c("h_z2(0)", "h_z2(1)")
    for (v in seq_along(folds)) {
        # train <- stack_data(origami::training(data, folds[[v]]), npsem)
        train <- stack_data(origami::training(data, folds[[v]]), varnames$M2)
        valid_1 <- valid_0 <- origami::validation(data, folds[[v]])
        valid_1[[varnames$A]] <- 1
        valid_0[[varnames$A]] <- 0
        
        p_zm2w <- crossfit(train, list(valid_0, valid_1), "tmp_tmce_delta",
                           c(varnames$Z, varnames$M2, varnames$W, varnames$A),
                           id = "tmp_tmce_id",
                           "binomial", learners = learners, bound = TRUE)
        
        p_m2w <- crossfit(train, list(valid_0, valid_1), "tmp_tmce_delta",
                          c(varnames$M2, varnames$W, varnames$A),
                          id = "tmp_tmce_id",
                          "binomial", learners = learners, bound = TRUE)
        
        h_z2[folds[[v]]$validation_set, "h_z2(0)"] <-
            (p_zm2w[[1]] / (1 - p_zm2w[[1]])) * ((1 - p_m2w[[1]]) / p_m2w[[1]])
        
        h_z2[folds[[v]]$validation_set, "h_z2(1)"] <-
            (p_zm2w[[2]] / (1 - p_zm2w[[2]])) * ((1 - p_m2w[[2]]) / p_m2w[[2]])
    }
    h_z2
}

h_z <- function(data, varnames, folds, learners, ...) {
    h_z <- matrix(nrow = nrow(data), ncol = 2)
    colnames(h_z) <- c("h_z(0)", "h_z(1)")
    for (v in seq_along(folds)) {
        # train <- stack_data(origami::training(data, folds[[v]]), npsem)
        train <- stack_data(origami::training(data, folds[[v]]), c(varnames$M1, varnames$M2))
        valid_1 <- valid_0 <- origami::validation(data, folds[[v]])
        valid_1[[varnames$A]] <- 1
        valid_0[[varnames$A]] <- 0
        
        p_zm1m2w <- crossfit(train, list(valid_0, valid_1), "tmp_tmce_delta",
                           c(varnames$Z, varnames$M1, varnames$M2, varnames$W, varnames$A),
                           id = "tmp_tmce_id",
                           "binomial", learners = learners, bound = TRUE)
        
        p_m1m2w <- crossfit(train, list(valid_0, valid_1), "tmp_tmce_delta",
                          c(varnames$M1, varnames$M2, varnames$W, varnames$A),
                          id = "tmp_tmce_id",
                          "binomial", learners = learners, bound = TRUE)
        
        h_z[folds[[v]]$validation_set, "h_z(0)"] <-
            (p_zm1m2w[[1]] / (1 - p_zm1m2w[[1]])) * ((1 - p_m1m2w[[1]]) / p_m1m2w[[1]])
        
        h_z[folds[[v]]$validation_set, "h_z(1)"] <-
            (p_zm1m2w[[2]] / (1 - p_zm1m2w[[2]])) * ((1 - p_m1m2w[[2]]) / p_m1m2w[[2]])
    }
    h_z
}


uniformly_sample_M <- function(data, M) {
    out <- foreach(m = M,
                   .combine = cbind,
                   .options.future = list(seed = TRUE)) %dofuture% {
                       sample_M(data[[m]])
                   }
    if (length(M) == 1) out <- as.matrix(out)
    colnames(out) <- M
    as.data.frame(out)
}




sample_M <- function(M) {
    vals_M <- unique(M)
    vals_M[sample.int(length(vals_M), size = length(M), replace = TRUE)]
}


stack_data <- function(data, Mnames) {
    delta <- data
    # delta[, npsem$M] <- uniformly_sample_M(data, npsem$M)[, npsem$M]
    delta[, Mnames] <- uniformly_sample_M(data, Mnames)[, Mnames]
    out <- rbind(data, delta)
    out[["tmp_tmce_delta"]] <- rep(c(0, 1), each = nrow(data))
    out[["tmp_tmce_id"]] <- rep(1:nrow(data), times = 2)
    out
}
