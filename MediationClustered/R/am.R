# models

# Treatment --------------------
a.c <- function(data_in, varnames, cluster_opt = "FE.glm", folds, learners, bounded = TRUE) {
    a_c <- matrix(nrow = nrow(data_in), ncol = 2)
    colnames(a_c) <- c("a(0|c)", "a(1|c)")
    v <- 1
    for (v in seq_along(folds)) {
        train <- origami::training(data_in, folds[[v]])
        valid <- origami::validation(data_in, folds[[v]])

        alist <- crossfit(train, list(valid), varnames$A, c(varnames$X), varnames,
                          ipw = NULL,
                          cluster_opt,
                          type = c("binomial"), learners, bounded)
        preds <- alist$preds
        a_c[folds[[v]]$validation_set, "a(0|c)"] <- 1 - preds[, 1]
        a_c[folds[[v]]$validation_set, "a(1|c)"] <- preds[, 1]
    }
    a_c
}


ipwA.c <- function(a_c, data_in, varnames, stablize = "none") {
  ipwA_c <- matrix(nrow = nrow(data_in), ncol = 2)
  colnames(ipwA_c) <- c("ipwA(0|c)", "ipwA(1|c)")
  ipwA_c[, "ipwA(0|c)"] <- (data_in[[varnames$A]]==0) * (1 / a_c[, "a(0|c)"])
  ipwA_c[, "ipwA(1|c)"] <- (data_in[[varnames$A]]==1) * (1 / a_c[, "a(1|c)"])
  if (stablize == "none") {
    sipwA_c <- ipwA_c
  }
  if (stablize == "non_cluster") {
    sipwA_c <- apply(ipwA_c, 2, function(w) { w / mean(w) })
  }
  if (stablize == "clustered") {
    clwA_c <- as.data.frame(ipwA_c) %>%
      mutate(S = data_in[[varnames$S]]) %>%
      group_by(S) %>%
      mutate(across(everything(), function(w) {
        if (mean(w) !=0 ) { clw <- w / mean(w)  }
        if (mean(w) ==0 ) { clw <- NA  }
        clw
        } ))
    sipwA_c <- as.matrix(clwA_c[, c("ipwA(0|c)", "ipwA(1|c)")])
  }

  sipwA_c
}



# Mediators --------------
# for one mediator p(M1=1|a,C)
m1.ac <- function(data_in, whichM = 1, varnames, ipw = NULL, cluster_opt = "FE.glm", folds, learners, bounded = FALSE) {
    m1_ac <- matrix(nrow = nrow(data_in), ncol = 2*2)
    colnames(m1_ac) <- c(glue("m{whichM}(0|0,c)"), glue("m{whichM}(1|0,c)"),
                         glue("m{whichM}(0|1,c)"), glue("m{whichM}(1|1,c)"))

    v <- 1
    for (v in seq_along(folds)) {
        train <- origami::training(data_in, folds[[v]])
        valid_1 <- valid_0 <- origami::validation(data_in, folds[[v]])
        valid_1[[varnames$A]] <- 1
        valid_0[[varnames$A]] <- 0

        if (length(ipw) > 0) {
          preds <- sapply(c(0, 1), function(a_val = 0) {
            valid_a <- origami::validation(data_in, folds[[v]])
            valid_a[[varnames$A]] <- a_val
            ind_a <- train[[varnames$A]] == a_val
            train_a <- train[ind_a, ]
            ipw_a <- ipw[ind_a, grep(a_val, colnames(ipw))]
            alist_a <- crossfit(train_a, list(valid_a),
                                varnames$M[whichM],
                                c(varnames$X),
                                varnames,
                                ipw_a,
                                cluster_opt,
                                type = c("binomial"), learners, bounded)
            alist_a$preds

          }, simplify = TRUE)

        }

        if (length(ipw) == 0) {
          alist <- crossfit(train, list(valid_0, valid_1),
                            varnames$M[whichM],
                            c(varnames$A, varnames$X),
                            varnames,
                            NULL,
                            cluster_opt,
                            type = c("binomial"), learners, bounded = FALSE)
          preds <- alist$preds
        }

        # valid_0
        m1_ac[folds[[v]]$validation_set, glue("m{whichM}(0|0,c)")] <- 1 - preds[, 1]
        m1_ac[folds[[v]]$validation_set, glue("m{whichM}(1|0,c)")] <- preds[, 1]
        # valid_1
        m1_ac[folds[[v]]$validation_set, glue("m{whichM}(0|1,c)")] <- 1 - preds[, 2]
        m1_ac[folds[[v]]$validation_set, glue("m{whichM}(1|1,c)")] <- preds[, 2]
    }
    m1_ac
}


m2.m1ac <- function(data_in, whichM = 2, varnames, ipw = NULL, cluster_opt = "FE.glm", interaction = c("AM.M1"), folds, learners, bounded = FALSE) {
    m2_m1ac <- matrix(nrow = nrow(data_in), ncol = 2+4*2)
    colnames(m2_m1ac) <- c(
      glue("m{whichM}(0|0,0,c)"), glue("m{whichM}(1|0,0,c)"),
                         glue("m{whichM}(0|1,0,c)"), glue("m{whichM}(1|1,0,c)"),
                         glue("m{whichM}(0|0,1,c)"), glue("m{whichM}(1|0,1,c)"),
                         glue("m{whichM}(0|1,1,c)"), glue("m{whichM}(1|1,1,c)"),
      glue("m{whichM}(0|obs,obs,c)"), glue("m{whichM}(1|obs,obs,c)")
      )

    v <- 1
    for (v in seq_along(folds)) {
        train <- origami::training(data_in, folds[[v]])
        valid <- valid_obs <- origami::validation(data_in, folds[[v]])
        valid_00 <- valid_10 <- valid_01 <- valid_11 <- origami::validation(data_in, folds[[v]])
        valid_00[, c(varnames$M[-whichM], varnames$A)] <- 0
        valid_10[, c(varnames$M[-whichM], varnames$A)] <- sapply(c(1, 0), rep, length = nrow(valid))
        valid_01[, c(varnames$M[-whichM], varnames$A)] <- sapply(c(0, 1), rep, length = nrow(valid))
        valid_11[, c(varnames$M[-whichM], varnames$A)] <- 1

        alist <- crossfit(train, list(valid_00, valid_10, valid_01, valid_11, valid_obs),
                          varnames$M[whichM],
                          c(varnames$M[-whichM], varnames$A, interaction, varnames$X),
                          varnames,
                          ipw,
                          cluster_opt,
                          type = c("binomial"), learners, bounded)
        preds <- alist$preds
        # valid_00
        m2_m1ac[folds[[v]]$validation_set, glue("m{whichM}(0|0,0,c)")] <- 1 - preds[, 1]
        m2_m1ac[folds[[v]]$validation_set, glue("m{whichM}(1|0,0,c)")] <- preds[, 1]
        # valid_10
        m2_m1ac[folds[[v]]$validation_set, glue("m{whichM}(0|1,0,c)")] <- 1 - preds[, 2]
        m2_m1ac[folds[[v]]$validation_set, glue("m{whichM}(1|1,0,c)")] <- preds[, 2]
        # valid_01
        m2_m1ac[folds[[v]]$validation_set, glue("m{whichM}(0|0,1,c)")] <- 1 - preds[, 3]
        m2_m1ac[folds[[v]]$validation_set, glue("m{whichM}(1|0,1,c)")] <- preds[, 3]
        # valid_11
        m2_m1ac[folds[[v]]$validation_set, glue("m{whichM}(0|1,1,c)")] <- 1 - preds[, 4]
        m2_m1ac[folds[[v]]$validation_set, glue("m{whichM}(1|1,1,c)")] <- preds[, 4]
        # valid_obs
        m2_m1ac[folds[[v]]$validation_set, glue("m{whichM}(0|obs,obs,c)")] <- 1 - preds[, 5]
        m2_m1ac[folds[[v]]$validation_set, glue("m{whichM}(1|obs,obs,c)")] <- preds[, 5]
    }
    m2_m1ac
}

m1m2.ac <- function(m1_ac, m2_m1ac) {
  m1m2_ac <- matrix(nrow = nrow(m1_ac), ncol = 8)
  vals <- expand.grid(m1 = c(0,1), m2 = c(0,1), a = c(0,1))
  colnames(m1m2_ac) <- glue("m1m2_ac({vals$m1},{vals$m2}|{vals$a},c)")

  for(i in 1:nrow(vals)) {
    m1 <- vals$m1[i]
    m2 <- vals$m2[i]
    a <- vals$a[i]
    m1m2_ac[, glue("m1m2_ac({m1},{m2}|{a},c)")] <-
      m2_m1ac[, glue("m2({m2}|{m1},{a},c)")] * m1_ac[, glue("m1({m1}|{a},c)")]
  }

  m1m2_ac
}

m2.ac <- function(m1m2_ac) {
  m2_ac <- matrix(nrow = nrow(m1m2_ac), ncol = 4)
  vals <- expand.grid(m2 = c(0,1), a = c(0,1))
  colnames(m2_ac) <- glue("m2({vals$m2}|{vals$a},c)")

  m1_vec <- c(0, 1)
  for (i in 1:nrow(vals)) {
    m2 <- vals$m2[i]
    a <- vals$a[i]
    m2_ac[, glue("m2({m2}|{a},c)")] <-
      rowSums( m1m2_ac[, glue("m1m2_ac({m1_vec},{m2}|{a},c)")] )
  }
  m2_ac
}

