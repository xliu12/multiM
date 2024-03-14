#' multiple mediators clustered data
#'
#' @export
#'
clustered <- function(data,
                      Sname,
                      Wnames = NULL, Xnames,
                      Aname,
                      Mnames, Mfamily = "binomial",
                      Yname, Yfamily = "binomial",
                      cluster_opt_a = "FE.glm",
                      cluster_opt_m = "FE.glm",
                      cluster_opt_y = "FE.glm",
                      interaction_fitm2 = c("AM.M1"),
                      interaction_fity = c("AM.M1", "AM.M2", "Mint", "AMint"),
                      num_folds = 1,
                      # bounds ,
                      learners = c("SL.glm") ,
                      learners_a = c("SL.glm"),
                      learners_m = c("SL.glm"),
                      learners_y = c("SL.glm")
) {
    # add interactions
    data_in <- data
    AM <- data_in[[Aname]] * data_in[, Mnames]
    colnames(AM) <- c("AM.M1", "AM.M2")
    Mint <- data.frame(data_in[, Mnames[1]] * data_in[, Mnames[2]])
    colnames(Mint) <- "Mint"
    AMint <- data.frame(data_in[[Aname]] * data_in[, Mnames[1]] * data_in[, Mnames[2]])
    colnames(AMint) <- "AMint"
    data_in <- data.frame(data, AM, Mint, AMint)
    # data_in <- data_in %>% mutate(
    #     AM = data_in[[Aname]] * data_in[, Mnames],
    #     Mint = data_in[, Mnames[1]] * data_in[, Mnames[2]],
    #     AMint = data_in[[Aname]] * data_in[, Mnames[1]] * data_in[, Mnames[2]]
    # )
    # data_in <- do.call(data.frame, data_in)

    Sdumm <- dummy_cols(data[[Sname]], remove_first_dummy = TRUE, remove_selected_columns = TRUE)
    colnames(Sdumm) <- paste0("S", 1:ncol(Sdumm))
    Sname_dummies <- colnames(Sdumm)

    data_in <- data.frame(data_in, Sdumm)
    data_in[[Sname]] <- as.factor(data_in[[Sname]])

    varnames <- list("A" = Aname, "M" = Mnames, "Y" = Yname,
                     "AM" = paste0(Aname, "M.", Mnames), "Mint" = "Mint", "AMint" = "AMint",
                     "S" = Sname, "Sdumm" = Sname_dummies,
                     "X" = Xnames, "W" = Wnames)

    Cnames <- list("X" = Xnames, "W" = Wnames, "Sdumm" = Sname_dummies)

    # num_folds <- 1 # no cross-fitting currently
    folds <- make_folds(cluster_ids = NULL, # Clusters are treated as a unit – that is, all observations within a cluster are placed in either the training or validation set.
                        strata_ids = data[[Sname]], # insofar as possible the distribution in the sample should be the same as the distribution in the training and validation sets.
                        fold_fun = folds_vfold, V = num_folds)
   

    if (num_folds==1) { folds[[1]]$training_set <- folds[[1]]$validation_set }
    # try
    # train <- valid <- data_in

    a_c <- a.c(data_in, varnames, cluster_opt_a, folds, learners_a, bounded = TRUE)
    # quantile(a_c)
    ipwA_c <- ipwA.c(a_c, data_in, varnames, stablize = "none") # later stablize in EIF weight

    if (Mfamily == "binomial") {
        m1_ac <- m1.ac(data_in, whichM = 1, varnames, ipw = NULL,
                       cluster_opt = cluster_opt_m, folds,
                       learners_m,
                       bounded = TRUE)

        m2_m1ac <- m2.m1ac(data_in, whichM = 2, varnames, ipw = NULL,
                           cluster_opt = cluster_opt_m,
                           interaction = interaction_fitm2,
                           folds, learners_m, bounded = TRUE)

        # joint probability p(M1,M2 | a, C)
        m1m2_ac <- m1m2.ac(m1_ac, m2_m1ac)
        m2_ac <- m2.ac(m1m2_ac)
    }

    y_m1m2ac <- y.m1m2ac(data_in, varnames, Yfamily = Yfamily, ipw = NULL, cluster_opt = cluster_opt_y, interaction = interaction_fity, folds, learners_y, bounded = FALSE)


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

    A <- data_in[[varnames$A]]
    M1 <- data_in[[varnames$M[1]]]
    M2 <- data_in[[varnames$M[2]]]
    Y <- data_in[[varnames$Y]]

    # multiply-robust
    thetas <- eifs <- list()
    # regression (G-computation)
    regs <- list()
    # weighting
    rmpw <- list()

    for (i in 1:nrow(a_vals)) {
        a0 <- a_vals$a0[i]
        a1 <- a_vals$a1[i]
        a2 <- a_vals$a2[i]
        ajo <- a_vals$ajo[i]

        # observed treatment
        ipwA <- ipwA_c[, glue("ipwA({a0}|c)")]
        # observed outcome
        mu <- y_m1m2ac[, "y_m1m2a(obs,obs,obs|c)"]
        # observed M distribution
        p_Mjo <- (M1==1)*(M2==1)*m1m2_ac[, glue("m1m2_ac(1,1|{a0},c)")] +
            (M1==1)*(M2==0)*m1m2_ac[, glue("m1m2_ac(1,0|{a0},c)")] +
            (M1==0)*(M2==1)*m1m2_ac[, glue("m1m2_ac(0,1|{a0},c)")] +
            (M1==0)*(M2==0)*m1m2_ac[, glue("m1m2_ac(0,0|{a0},c)")]
        # quantile(p_Mjo, c(0, 0.1, 0.9))
        # den <- a_c[, glue("a({a0}|c)")] * p_Mjo
        # quantile(den)

        if (!is.na(a1)) {
            p_M1 <- (M1==1)*m1_ac[, glue("m1(1|{a1},c)")] + (M1==0)*m1_ac[, glue("m1(0|{a1},c)")]
            p_M2 <- (M2==1)*m2_ac[, glue("m2(1|{a2},c)")] + (M2==0)*m2_ac[, glue("m2(0|{a2},c)")]

            # integrate over p(m2 | a2, c)
            mu_M1a0c <- y.M1a0c(y_m1m2ac, a0, m2_ac, a2)
            # integrate over p(m1 | a1, c)
            mu_M2a0c <- y.M2a0c(y_m1m2ac, a0, m1_ac, a1)
            # integrate over p(m1 | a1, c)*p(m2 | a2, c)
            mu_a0c_12 <- y.a0c_12(y_m1m2ac, a0, m1_ac, a1, m2_ac, a2)

            # fot the effect via M1 or M2 alone ----
            h_M1M2 <- p_M1 * p_M2 / p_Mjo
            eify <- 1*(A == a0)*ipwA*h_M1M2 / mean(1*(A == a0)*ipwA*h_M1M2) * (Y - mu)

            ipwa1 <- ipwA_c[, glue("ipwA({a1}|c)")]
            eifm1 <- 1*(A == a1)*ipwa1 / mean(1*(A == a1)*ipwa1) * (mu_M1a0c - mu_a0c_12)

            ipwa2 <- ipwA_c[, glue("ipwA({a2}|c)")]
            eifm2 <- 1*(A == a2)*ipwa2 / mean(1*(A == a2)*ipwa2) * (mu_M2a0c - mu_a0c_12)

            # eif with marginal mediator distributions
            eif_mar <- eify + eifm1 + eifm2 + mu_a0c_12

            eifs[[glue("Y({a0},gm1({a1}),gm2({a2}))")]] <- eif_mar
            thetas[[glue("Y({a0},gm1({a1}),gm2({a2}))")]] <- mean(eif_mar)

            rmpw[[glue("Y({a0},gm1({a1}),gm2({a2}))")]] <- mean( 1*(A == a0)*ipwA*h_M1M2 / mean(1*(A == a0)*ipwA*h_M1M2) * Y )
            regs[[glue("Y({a0},gm1({a1}),gm2({a2}))")]] <- mean(mu_a0c_12)

        }

        # for the joint ajo ----
        if (!is.na(ajo)) {
            p_Majo <- (M1==1)*(M2==1)*m1m2_ac[, glue("m1m2_ac(1,1|{ajo},c)")] +
                (M1==1)*(M2==0)*m1m2_ac[, glue("m1m2_ac(1,0|{ajo},c)")] +
                (M1==0)*(M2==1)*m1m2_ac[, glue("m1m2_ac(0,1|{ajo},c)")] +
                (M1==0)*(M2==0)*m1m2_ac[, glue("m1m2_ac(0,0|{ajo},c)")]
            # integrate over p(m1, m2 | ajo, c)
            mu_a0c_jo <- y.a0c_jo(y_m1m2ac, a0, m1m2_ac, ajo)

            h_Mjo <- p_Majo / p_Mjo
            # (1/ipwA)*p_Mjo %>% quantile(., c(0, 0.1, 0.9, 0.99))
            eify_ajo <- 1*(A == a0)*ipwA*h_Mjo / mean(1*(A == a0)*ipwA*h_Mjo) * (Y - mu)
            # fixing A = a0
            mu_M1M2a0c <- y_m1m2ac[, glue("y_m1m2a(obs,obs,{a0}|c)")]
            ipwajo <- ipwA_c[, glue("ipwA({ajo}|c)")]
            eifm_ajo <- 1*(A == ajo)*ipwajo / mean(1*(A == ajo)*ipwajo) * (mu_M1M2a0c - mu_a0c_jo)

            # eif with the joint mediator distribution
            eif_jo <- eify_ajo + eifm_ajo + mu_a0c_jo

            eifs[[glue("Y({a0},gmjo({ajo}))")]] <- eif_jo
            thetas[[glue("Y({a0},gmjo({ajo}))")]] <- mean(eif_jo)

            rmpw[[glue("Y({a0},gmjo({ajo}))")]] <- mean( 1*(A == a0)*ipwA*h_Mjo / mean(1*(A == a0)*ipwA*h_Mjo) * (Y) )
            regs[[glue("Y({a0},gmjo({ajo}))")]] <- mean(mu_a0c_jo)
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

    # regression
    regs_effs <- unlist(effect(regs))
    # weighting
    rmpw_effs <- unlist(effect(rmpw))

    out <- mget(ls(envir = environment()))

    return(out)
}



effect <- function(thetas) {
    thetas[["IDE(,gmjo(0))"]] <- thetas[["Y(1,gmjo(0))"]] - thetas[["Y(0,gmjo(0))"]]
    thetas[["IIE_Mjo(1,)"]] <- thetas[["Y(1,gmjo(1))"]] - thetas[["Y(1,gmjo(0))"]]
    thetas[["IIE_M1(1,,gm2(0))"]] <- thetas[["Y(1,gm1(1),gm2(0))"]] - thetas[["Y(1,gm1(0),gm2(0))"]]
    thetas[["IIE_M2(1,gm1(1),)"]] <- thetas[["Y(1,gm1(1),gm2(1))"]] - thetas[["Y(1,gm1(1),gm2(0))"]]
    thetas[["IIE_Mdep"]] <- thetas[["IIE_Mjo(1,)"]] - thetas[["IIE_M1(1,,gm2(0))"]] - thetas[["IIE_M2(1,gm1(1),)"]]

    thetas
}
