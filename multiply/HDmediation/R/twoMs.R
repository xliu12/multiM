not_transported <- function(data, Aname, Wnames, Znames, M1names, M2names, Yname, family, num_folds = 1, partial_tmle, bounds,
                            learners_g = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                            learners_e = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                            learners_b = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                            learners_hz = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                            learners_u = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                            learners_ubar = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                            learners_v = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                            learners_vbar = c("SL.glm", "SL.glm.interaction", "SL.mean")) {
    # npsem <- Npsem$new(A = A, W = W, Z = Z, M = M, Y = Y)
    varnames <- list(A = Aname, W = Wnames, Z = Znames, M1 = M1names, M2 = M2names, Y = Yname)
    # varnames <- list(tt = Trtname, C = Cnames, Z = Znames, M1 = M1names, M2 = M2names, Y = Yname)
    
    folds <- make_folds(data, num_folds)

    bounds <- scale_y(data[[varnames$Y]], family, bounds)
    data[[varnames$Y]] <- Y <- bounds$y
    A <- data[[varnames$A]]

    gg <- p_t.c(data, varnames, folds, learners_g)
    em1 <- p_t.m1c(data, varnames, folds, learners_e)
    em2 <- p_t.m2c(data, varnames, folds, learners_e)
    em1m2 <- p_t.m1m2c(data, varnames, folds, learners_e)
    bb <- mu(data, varnames, family, folds, learners_b)
    
    vvbar <- vector("list", length = 6)
    names(vvbar) <- c("Y(0,M(t_jo=0))", "Y(1,M(t_jo=0))", "Y(1,M(t_jo=1))", 
                      "Y(1,M1(1),M2(1))", "Y(1,M1(0),M2(1))", "Y(1,M1(0),M2(0))")
    eifs <- list()
    thetas <- list()
    
    for (param in list(c(t0 = 0, t_jo = 0, t_1 = NA, t_2 = NA), 
                       c(t0 = 1, t_jo = 0, t_1 = NA, t_2 = NA), 
                       c(t0 = 1, t_jo = 1, t_1 = NA, t_2 = NA),
                       c(t0 = 1, t_jo = NA, t_1 = 1, t_2 = 1),
                       c(t0 = 1, t_jo = NA, t_1 = 0, t_2 = 1),
                       c(t0 = 1, t_jo = NA, t_1 = 0, t_2 = 0)
                       )) {
        
        t0 <- param["t0"]
        
        if (!is.na(param["t_jo"])) {
            t_jo <- param["t_jo"]
            hz <- h_z(data, varnames, folds, learners_hz)
            hmjo <- h_mjo(hz, gg, em1m2, t_jo, t0)  
            muMjo <- mu_Mjo(data, varnames, bb, hmjo, t0, folds, learners_u)
            mubarMjo <- mubar_Mjo(data, varnames, muMjo, t0, folds, learners_ubar)
            muZ <- mu_Z(data, varnames, bb, hz, t0, folds, learners_v)
            mubarZ <- mubar_Z(data, varnames, muZ, t_jo, folds, learners_vbar)
            
            # EIF calculation
            ipwy <- (A == t0) / gg[, gl("g({t0}|w)")]
            eify <- ipwy * hmjo / mean(ipwy * hmjo) * (Y - bb[, gl("mu({t0},Z,M1,M2,C)")])
            
            ipwz <- (A == t0) / gg[, gl("g({t0}|w)")]
            eifz <- ipwz / mean(ipwz) * (muMjo[, 1] - mubarMjo[, 1])
            
            ipwm <- (A == t_jo) / gg[, gl("g({t_jo}|w)")]
            eifm <- ipwm / mean(ipwm) * (muZ[, 1] - mubarZ[, 1])
            
            eif <- rescale_y(eify + eifz + eifm + mubarZ[, 1], bounds$bounds)
            theta <- mean(eif)
            
            vname <- paste0("Y(", t0 ,",M(t_jo=", t_jo, "))")
            vvbar[[vname]] <- mubarZ[, 1]
            
            thetas[vname] <- theta
            eifs[[vname]] <- eif
        }
        
        if (is.na(param[["t_jo"]])) {
            t_1 <- param["t_1"]
            t_2 <- param["t_2"]
            hz1 <- h_z1(data, varnames, folds, learners_hz)
            hz2 <- h_z2(data, varnames, folds, learners_hz)
            hzm1 <- h_zm1(data, varnames, folds, learners_hz)
            hzm2 <- h_zm2(data, varnames, folds, learners_hz)
            
            hm1 <- h_m1(hz1, gg, em1, t_1, t0)
            hm2sec <- h_m2sec(hzm1, gg, em2, t_2, t0)
            hm1sec <- h_m1sec(hzm2, gg, em1, t_1, t0)
            
            hm1m2 <- hm1 * hm2sec
            hm1z <- hm1sec * hz2[, gl("h_z2({t0})")] 
            hm2z <- hm2sec * hz1[, gl("h_z1({t0})")]
            
            muM1M2 <- mu_M1M2(data, varnames, bb, hm1m2, t0, folds, learners_u)
            mubarM1M2 <- mubar_M1M2(data, varnames, muM1M2, t0, folds, learners_ubar)
            
            muZM1 <- mu_ZM1(data, varnames, bb, hm1z, t0, folds, learners_v)
            mubarZM1 <- mubar_ZM1(data, varnames, muZM1, t_2, folds, learners_vbar)
            
            muZM2 <- mu_ZM2(data, varnames, bb, hm2z, t0, folds, learners_v)
            mubarZM2 <- mubar_ZM2(data, varnames, muZM2, t_1, folds, learners_vbar)
            
            # EIF calculation
            ipwy <- (A == t0) / gg[, gl("g({t0}|w)")]
            eify <- ipwy * hm1m2 / mean(ipwy * hm1m2) * (Y - bb[, gl("mu({t0},Z,M1,M2,C)")])
            
            ipwz <- (A == t0) / gg[, gl("g({t0}|w)")]
            eifz <- ipwz / mean(ipwz) * (muM1M2[, 1] - mubarM1M2[, 1])
            
            ipwm1 <- (A == t_1) / gg[, gl("g({t_1}|w)")]
            eifm1 <- ipwm1 / mean(ipwm1) * (muZM2[, 1] - mubarZM2[, 1])
            
            ipwm2 <- (A == t_2) / gg[, gl("g({t_2}|w)")]
            eifm2 <- ipwm2 / mean(ipwm2) * (muZM1[, 1] - mubarZM1[, 1])
            
            eif <- rescale_y(eify + eifz + eifm1 + eifm2 + mubarZM1[, 1], bounds$bounds)
            theta <- mean(eif)
            
            vname <- paste0("Y(", t0 ,",M1(", t_1, "),", "M2(", t_2, ")", ")")
            vvbar[[vname]] <- mubarZM1[, 1]
            
            thetas[vname] <- theta
            eifs[[vname]] <- eif
        }
        
        partial_tmle <- FALSE
        if (partial_tmle) {
            fit <- glm(Y ~ 1, offset = qlogis(bb[, gl("b({aprime},Z,M,W)")]), family = "binomial",
                       subset = A == aprime, weights = ipwy * hm / mean(ipwy * hm))
            bb[, gl("b({aprime},Z,M,W)")] <- plogis(coef(fit) + qlogis(bb[, gl("b({aprime},Z,M,W)")]))
        }

    }

    
    
    
    # multiply-robust 
    eifs_effs <- effect(eifs)
    thetas_effs <- sapply(eifs_effs, mean)
    se_effs <- sapply(eifs_effs, function(s) {
        sqrt(var(s) / nrow(data))
    })
    ci_effs <- sapply(eifs_effs, function(s) {
        mean(s) + c(-1, 1) * qnorm(0.975) * sqrt(var(s) / nrow(data))
    })
    
    # sequential regression
    vvbar_effs <- effect(vvbar)
    reg_effs <- sapply(vvbar_effs, mean, na.rm = TRUE)
    
    # output ------------
    estimates <- data.frame(names(thetas_effs), 
               cbind(thetas_effs, se_effs, t(ci_effs), reg_effs))
    colnames(estimates) <- c("effect", "est_onestep", "se", "ci1", "ci2", "est_reg")
    
    # out <- list(estimates = estimates, eifs_effs = eifs_effs)
    estimates
}



effect <- function(thetas) {
    thetas[["IDE(,M(t_jo=0))"]] <- thetas[["Y(1,M(t_jo=0))"]] - thetas[["Y(0,M(t_jo=0))"]]
    thetas[["IIE_Mjo(1,)"]] <- thetas[["Y(1,M(t_jo=1))"]] - thetas[["Y(1,M(t_jo=0))"]]
    thetas[["IIE_M1(1,,M2(1))"]] <- thetas[["Y(1,M1(1),M2(1))"]] - thetas[["Y(1,M1(0),M2(1))"]]
    thetas[["IIE_M2(1,M1(0),)"]] <- thetas[["Y(1,M1(0),M2(1))"]] - thetas[["Y(1,M1(0),M2(0))"]]
    thetas[["IIE_Mdep"]] <- thetas[["Y(1,M(t_jo=1))"]] - thetas[["Y(1,M1(1),M2(1))"]] - (thetas[["Y(1,M(t_jo=0))"]] - thetas[["Y(1,M1(0),M2(0))"]])
    
    thetas
}
