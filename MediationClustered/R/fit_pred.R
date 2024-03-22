# fit a model


crossfit <- function(train, valid.list, yname, xnames, varnames,
                     ipw = NULL,
                     cluster_opt = "FE.glm",
                     type, learners, bounded = FALSE) {

    Sname <- varnames$S
    Sname_dummies <- varnames$Sdumm

    family <- ifelse(type == "binomial", binomial(), gaussian())

    df_lm <- data.frame(Y = train[[yname]],
                        train[, c(xnames, varnames$W), drop = FALSE],
                        S = train[[Sname]] )
    if (length(ipw) > 0) {
        df_lm$wreg <- ipw
    }
    if (length(ipw) == 0) {
        df_lm$wreg <- rep(1, nrow(df_lm))
    }
    # fixed effects
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

    if (cluster_opt %in% c("fix.mlr", "noncluster.mlr")) {
        df_FE <- data.frame(df_lm, train[, Sname_dummies])
        set.seed(12345)
        if (cluster_opt == "fix.mlr") {

            fit <- SuperLearner::SuperLearner(
                df_FE$Y,
                df_FE[, c(xnames, varnames$W, Sname_dummies), drop = FALSE],
                obsWeights = df_FE$wreg,
                family = family[[1]],
                SL.library = learners
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
            newX <- data.frame(validX[, c(xnames, varnames$W, Sname_dummies), drop = FALSE])
            preds <- predict(fit, newX[, fit$varNames])$pred
            if (!bounded) {
                return(preds)
            }
            bound(preds)
        }, simplify = TRUE)
    }

    # random effects
    if (cluster_opt == "RE.glm") {
        REformula <- paste("Y ~", paste(c(xnames, varnames$W), collapse = " + "), "+ (1 | S)")
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
            newX <- data.frame(validX[, c(xnames, varnames$W), drop = FALSE],
                               S = validX[[Sname]] )
            preds <- predict(fit, newX, type = "response")
            if (!bounded) {
                return(preds)
            }
            bound(preds)
        }, simplify = TRUE)

    }


    # clusters ------------------------
    if (cluster_opt == "FE.mlr") {
        df_FE <- data.frame(df_lm, train[, Sname_dummies])
        df_cwc <- df_lm %>%
            group_by(S) %>%
            mutate(across(everything(), ~.x-mean(.x) ))
        df_cluster <- data.frame(df_cwc, train[, Sname_dummies])
        fit <- SuperLearner::SuperLearner(
            df_FE$Y,
            df_cluster[, c(xnames, Sname_dummies), drop = FALSE],
            family = "gaussian",
            SL.library = learners,
            env = environment(SuperLearner::SuperLearner)
        )
        preds <- sapply(valid.list, function(validX) {
            newX <- data.frame(Y = validX[[yname]],
                               validX[, xnames, drop = FALSE],
                               S = validX[[Sname]] )
            newX_cwc <- newX %>%
                group_by(S) %>%
                mutate(across(everything(), ~.x-mean(.x) ))
            newX_cluster <- data.frame(newX_cwc, validX[, Sname_dummies])

            preds <- predict(fit, newX_cluster[, fit$varNames])$pred

            if (!bounded) {
                return(preds)
            }
            bound(preds)
        }, simplify = TRUE)

    }


    # demean
    if (cluster_opt == "demean") { # continuous Y
        df_FE <- data.frame(df_lm, train[, Sname_dummies])
        df_demean <- df_lm %>%
            group_by(S) %>%
            mutate(across(everything(), ~.x-mean(.x) ))
        df_demean_FE <- data.frame(df_demean, train[, Sname_dummies])
        fit <- SuperLearner::SuperLearner(
            df_demean_FE$Y,
            df_demean_FE[, c(xnames, Sname_dummies), drop = FALSE],
            family = "gaussian",
            SL.library = learners,
            env = environment(SuperLearner::SuperLearner)
        )
        preds <- sapply(valid.list, function(validX) {
            newX <- data.frame(Y = validX[[yname]],
                               validX[, xnames, drop = FALSE],
                               S = validX[[Sname]] )
            newX_demean <- newX %>%
                group_by(S) %>%
                mutate(across(everything(), ~.x-mean(.x) ))
            newX_demean_FE <- data.frame(newX_demean, validX[, Sname_dummies])
            preds_demean <- predict(fit, newX_demean_FE[, fit$varNames])$pred
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



bound <- function(vals, tol = 1e-6) {
    vals[vals < tol] <- tol
    vals[vals > 1 - tol] <- 1 - tol
    return(vals)
}



