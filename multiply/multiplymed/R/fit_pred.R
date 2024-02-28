# fit a nuisance model
library(lme4)
library(fastDummies)



crossfit <- function(train, valid.list, yname, xnames, varnames,
                     ipw = NULL, 
                     cluster_opt = "FE.glm",  
                     type, learners, bounded = FALSE) {
    
    Sname <- varnames$S
    Sname_dummies <- varnames$Sdumm
    # clust_id <- train[[Sname]]
    family <- ifelse(type == "binomial", binomial(), gaussian())
    
    df_lm <- data.frame(Y = train[[yname]], 
                        train[, xnames, drop = FALSE], 
                        S = train[[Sname]] )
    if (length(ipw) > 0) {
        df_lm$wreg <- ipw
    }
    if (length(ipw) == 0) {
        df_lm$wreg <- rep(1, nrow(df_lm))
    }
    # fixed effects -----------------
    if (cluster_opt %in% c("FE.glm", "noncluster.glm")) {
        if (cluster_opt == "FE.glm") {
            fit <- glm(formula = paste("Y ~ S +", paste(xnames, collapse = " + ")), 
                       weights = wreg, 
                       data = df_lm, family = family[[1]])
        }
        if (cluster_opt == "noncluster.glm") {
            fit <- glm(formula = paste("Y ~", paste(xnames, collapse = " + ")), 
                       weights = wreg, 
                       data = df_lm, family = family[[1]])
        }
        preds <- sapply(valid.list, function(validX) {
            newX <- data.frame(validX[, xnames, drop = FALSE], 
                               S = validX[[Sname]] )
            preds <- predict(fit, newX, type = "response")
            if (!bounded) {
                return(preds)
            }
            bound(preds)
        }, simplify = TRUE)
    }
    
    if (cluster_opt %in% c("FE.mlr", "noncluster.mlr")) {
        df_FE <- data.frame(df_lm, train[, Sname_dummies])
        
        if (cluster_opt == "FE.mlr") {
            fit <- SuperLearner::SuperLearner(
                df_FE$Y,
                df_FE[, c(xnames, Sname_dummies), drop = FALSE],
                obsWeights = df_FE$wreg,
                family = family[[1]],
                SL.library = learners
                # id = clust_id, # Optional cluster identification variable. For the cross-validation splits, id forces observations in the same cluster to be in the same validation fold.
            )
        }
        if (cluster_opt == "noncluster.mlr") {
            fit <- SuperLearner::SuperLearner(
                df_FE$Y,
                df_FE[, c(xnames), drop = FALSE],
                obsWeights = df_FE$wreg,
                family = family[[1]],
                SL.library = learners
            )
        }
        
        preds <- sapply(valid.list, function(validX) {
            newX <- data.frame(validX[, c(xnames, Sname_dummies), drop = FALSE])
            preds <- predict(fit, newX[, fit$varNames])$pred 
            if (!bounded) {
                return(preds)
            }
            bound(preds)
        }, simplify = TRUE)
    }
    
    # random effects ------------------------
    if (cluster_opt == "RE.glm") {
        REformula <- paste("Y ~", paste(xnames, collapse = " + "), "+ (1 | S)")
        if (family[[1]] == "gaussian") {
            fit <- lmer(formula = REformula, 
                         weights = wreg,
                         data = df_lm) 
        }
        if (family[[1]] != "gaussian") {
            fit <- glmer(formula = REformula, 
                        weights = wreg,
                        data = df_lm, family = family[[1]]) 
        }
        preds <- sapply(valid.list, function(validX) {
            newX <- data.frame(validX[, xnames, drop = FALSE], 
                               S = validX[[Sname]] )
            preds <- predict(fit, newX, type = "response")
            if (!bounded) {
                return(preds)
            }
            bound(preds)
        }, simplify = TRUE)
        
    }
    # demean ------------------------
    if (cluster_opt == "demean") { # continuous Y
        df_demean <- df_lm %>% 
            group_by(S) %>% 
            mutate(across(everything(), ~.x-mean(.x) ))
        # df_demean %>% group_by(S) %>% summarise(across(everything(), ~mean(.x) ))
        fit <- SuperLearner::SuperLearner(
            df_demean$Y,
            df_demean[, xnames, drop = FALSE],
            obsWeights = wreg,
            family = "gaussian",
            SL.library = learners,
            # id = clust_id, # Optional cluster identification variable. For the cross-validation splits, id forces observations in the same cluster to be in the same validation fold.
            env = environment(SuperLearner::SuperLearner)
        )
        preds <- sapply(valid.list, function(validX) {
            newX <- data.frame(Y = validX[[yname]], 
                validX[, xnames, drop = FALSE], 
                               S = validX[[Sname]] )
            newX_demean <- newX %>% 
                group_by(S) %>% 
                mutate(across(everything(), ~.x-mean(.x) ))
            
            preds_demean <- predict(fit, newX_demean[, fit$varNames])$pred
            # add back the cluster means of Y
            newXcm <- newX %>% 
                group_by(S) %>% 
                mutate(Y_cm = mean(Y)) 
            preds <- preds_demean + newXcm$Y_cm
            if (!bounded) {
                return(preds)
            }
            bound(preds)
        }, simplify = TRUE)
        
    }
    
    out <- list(fit = fit, preds = preds)
    return(out)
}

# yname <- "A"
# xnames <- Xnames


bound <- function(vals, tol = 1e-6) {
    vals[vals < tol] <- tol
    vals[vals > 1 - tol] <- 1 - tol
    return(vals)
}



# SL.modified ----------------

SL.glm.saturated <- function(Y, X, newX, family, obsWeights, ...) {
    if (is.matrix(X)) {
        X = as.data.frame(X)
    }
    f <- as.formula(paste0("Y ~ .^", ncol(X)))
    fit.glm <- glm(f, data = X, family = family, weights = obsWeights)
    if (is.matrix(newX)) {
        newX = as.data.frame(newX)
    }
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred = pred, fit = fit)
    return(out)
}

#' @export
SL.lightgbm <- function(Y, X, newX, family, obsWeights, id, nrounds = 1000, verbose = -1,
                        learning_rate = 0.1, min_data_in_leaf = 10, max_depth = -1, ...) {
    if (!requireNamespace("lightgbm", quietly = FALSE)) {
        stop("loading required package (lightgbm) failed", call. = FALSE)
    }
    
    if (family$family == "gaussian") {
        objective <- "regression"
        evalu <- ""
    }
    
    if (family$family == "binomial") {
        objective <- "binary"
        evalu <- "binary_logloss"
    }
    
    if (!is.matrix(X)) {
        X <- model.matrix(~. - 1, X)
    }
    
    lgb_data <- try(
        lightgbm::lgb.Dataset(
            data = X,
            label = as.numeric(Y)
        ), silent = TRUE
    )
    
    try(lightgbm::set_field(lgb_data, "weight", as.numeric(obsWeights)), silent = TRUE)
    
    params <- list(
        min_data_in_leaf = min_data_in_leaf,
        learning_rate = learning_rate,
        max_depth = max_depth
    )
    
    model <- lightgbm::lgb.train(params, data = lgb_data, obj = objective, eval = evalu,
                                 nrounds = nrounds, verbose = verbose)
    
    if (!is.matrix(newX)) {
        newX <- model.matrix(~. - 1, newX)
    }
    
    pred <- predict(model, newX)
    fit <- list(object = model)
    class(fit) <- c("SL.lightgbm")
    out <- list(pred = pred, fit = fit)
    return(out)
}

#' @export
predict.SL.lightgbm <- function(object, newdata, family, ...) {
    if (!requireNamespace("lightgbm", quietly = FALSE)) {
        stop("loading required package (lightgbm) failed", call. = FALSE)
    }
    
    if (!is.matrix(newdata)) {
        newdata <- model.matrix(~. - 1, newdata)
    }
    pred <- predict(object$object, newdata)
    return(pred)
}
