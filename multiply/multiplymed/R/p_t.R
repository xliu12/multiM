
p_t.c <- function(data, varnames, folds, learners) {
    gmat <- matrix(nrow = nrow(data), ncol = 2)
    colnames(gmat) <- c("g(0|w)", "g(1|w)")
    
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])

        preds <- crossfit(train, list(valid), varnames$A, c(varnames$W),
                          "binomial", learners = learners, bound = TRUE)[[1]]
        
        gmat[folds[[v]]$validation_set, "g(0|w)"] <- 1 - preds
        gmat[folds[[v]]$validation_set, "g(1|w)"] <- preds
    }
    gmat
}



p_t.m1c <- function(data, varnames, folds, learners) {
    emat <- matrix(nrow = nrow(data), ncol = 2)
    colnames(emat) <- c("e(0|m1,w)", "e(1|m1,w)")
    
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        
        preds <- crossfit(train, list(valid), varnames$A, c(varnames$M1, varnames$W),
                          "binomial", learners = learners, bound = TRUE)[[1]]
        emat[folds[[v]]$validation_set, "e(0|m1,w)"] <- 1 - preds
        emat[folds[[v]]$validation_set, "e(1|m1,w)"] <- preds
    }
    
    em1 <- emat
    em1
}


p_t.m2c <- function(data, varnames, folds, learners) {
    emat <- matrix(nrow = nrow(data), ncol = 2)
    colnames(emat) <- c("e(0|m2,w)", "e(1|m2,w)")
    
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        
        preds <- crossfit(train, list(valid), varnames$A, c(varnames$M2, varnames$W),
                          "binomial", learners = learners, bound = TRUE)[[1]]
        emat[folds[[v]]$validation_set, "e(0|m2,w)"] <- 1 - preds
        emat[folds[[v]]$validation_set, "e(1|m2,w)"] <- preds
    }
    
    em2 <- emat
    em2
}


p_t.m1m2c <- function(data, varnames, folds, learners) {
    emat <- matrix(nrow = nrow(data), ncol = 2)
    colnames(emat) <- c("e(0|m1,m2,w)", "e(1|m1,m2,w)")
    
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        
        preds <- crossfit(train, list(valid), varnames$A, c(varnames$M1, varnames$M2, varnames$W),
                          "binomial", learners = learners, bound = TRUE)[[1]]
        emat[folds[[v]]$validation_set, "e(0|m1,m2,w)"] <- 1 - preds
        emat[folds[[v]]$validation_set, "e(1|m1,m2,w)"] <- preds
    }
    
    em1m2 <- emat
    em1m2
}
