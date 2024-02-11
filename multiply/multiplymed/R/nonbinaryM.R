# (Not now) if m1 or m2 are not binary
if (Mfamily != "binomial") {
    a_m1c <- a.m1c(data_in, whichM = 1, varnames, cluster_opt = "FE", folds, learners)
    a_m2c <- a.m1c(data_in, whichM = 2, varnames, cluster_opt = "FE", folds, learners)
    a_mc <- a.mc(data_in, varnames, cluster_opt = "FE", folds, learners)
}


gg <- p_t.c(data_in, varnames, folds, learners_g)
em1 <- p_t.m1c(data_in, varnames, folds, learners_e)
em2 <- p_t.m2c(data_in, varnames, folds, learners_e)
em1m2 <- p_t.m1m2c(data_in, varnames, folds, learners_e)
bb <- mu(data_in, varnames, family, folds, learners_b)

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
        hz <- h_z(data_in, varnames, folds, learners_hz)
        hmjo <- h_mjo(hz, gg, em1m2, t_jo, t0)  
        muMjo <- mu_Mjo(data_in, varnames, bb, hmjo, t0, folds, learners_u)
        mubarMjo <- mubar_Mjo(data_in, varnames, muMjo, t0, folds, learners_ubar)
        muZ <- mu_Z(data_in, varnames, bb, hz, t0, folds, learners_v)
        mubarZ <- mubar_Z(data_in, varnames, muZ, t_jo, folds, learners_vbar)
        
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
        hz1 <- h_z1(data_in, varnames, folds, learners_hz)
        hz2 <- h_z2(data_in, varnames, folds, learners_hz)
        hzm1 <- h_zm1(data_in, varnames, folds, learners_hz)
        hzm2 <- h_zm2(data_in, varnames, folds, learners_hz)
        
        hm1 <- h_m1(hz1, gg, em1, t_1, t0)
        hm2sec <- h_m2sec(hzm1, gg, em2, t_2, t0)
        hm1sec <- h_m1sec(hzm2, gg, em1, t_1, t0)
        
        hm1m2 <- hm1 * hm2sec
        hm1z <- hm1sec * hz2[, gl("h_z2({t0})")] 
        hm2z <- hm2sec * hz1[, gl("h_z1({t0})")]
        
        muM1M2 <- mu_M1M2(data_in, varnames, bb, hm1m2, t0, folds, learners_u)
        mubarM1M2 <- mubar_M1M2(data_in, varnames, muM1M2, t0, folds, learners_ubar)
        
        muZM1 <- mu_ZM1(data_in, varnames, bb, hm1z, t0, folds, learners_v)
        mubarZM1 <- mubar_ZM1(data_in, varnames, muZM1, t_2, folds, learners_vbar)
        
        muZM2 <- mu_ZM2(data_in, varnames, bb, hm2z, t0, folds, learners_v)
        mubarZM2 <- mubar_ZM2(data_in, varnames, muZM2, t_1, folds, learners_vbar)
        
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


