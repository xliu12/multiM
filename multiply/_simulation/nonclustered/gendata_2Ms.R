w_n <- function(n = 200, num_w = 3) { # w: baseline covariates
  rmvnorm(n, mean = rep(0, num_w))
}
g <- function(a, w, randomized = FALSE) {
  pscore <- pnorm( rowSums(sqrt(0.13/ncol(w)) * w) )
  if (randomized) {pscore <- .5} 
  a * pscore + (1 - a) * (1 - pscore)
}

# pz1 <- function(z, a, w, binary = TRUE) {
#   aw <- 0.1*a + rowSums(sqrt(0.13/ncol(w)) * w)
#   if (!binary) {
#     aw + rnorm(nrow(w), sd = sqrt(1-0.13))
#   }
#   if (binary) { 
#     prob1 <- pnorm(aw)
#     z * prob1 + (1 - z) * (1 - prob1)
#   }
# }
# 
# pz2 <- function(z, a, w, binary = TRUE) {
#   aw <- 0.1*a + rowSums(sqrt(0.13/ncol(w)) * w)
#   if (!binary) {
#     aw + rnorm(nrow(w), sd = sqrt(1-0.13))
#   }
#   if (binary) { 
#     prob1 <- pnorm(aw)
#     z * prob1 + (1 - z) * (1 - prob1)
#   }
# }
# 
# pm1 <- function(m, z1, a, w) {
#     prob1 <- 0.6 + 0.1*z1 + 0.05*a - 0.3*w
#     m * prob1 + (1 - m) * (1 - prob1)
# }
# 
# pm2 <- function(m, z2, a, w) {
#     prob1 <- 0.33 + 0.22*z2 + 0.05*a + 0.15*w
#     m * prob1 + (1 - m) * (1 - prob1)
# }

pz <- function(z, a, w, binary = TRUE) {
  aw <- 0.1*a + rowSums(sqrt(0.13/ncol(w)) * w)
  if (!binary) {
    out <- aw + rnorm(nrow(w), sd = sqrt(1-0.13))
  }
  if (binary) { 
    prob1 <- pnorm(aw)
    out <- z * prob1 + (1 - z) * (1 - prob1)
  }
  out
}

pm <- function(m1, m2, z, a, w, binary = TRUE) {
  U <- rnorm(nrow(w)) # common causes of mediators
  Em1 <- 0.1*z + 0.1*a + 0.2*U + rowSums(sqrt(0.13/ncol(w)) * w)
  Em2 <- 0.3*z + 0.3*a + 0.2*U + rowSums(sqrt(0.13/ncol(w)) * w)
  
  if (!binary) {
    em <- rmvnorm(nrow(w), mean = c(0,0), sigma = diag((1-0.2^2-0.13), nrow = 2))
    out <- cbind(Em1, Em2) + em
  }
  if (binary) { 
    prob1 <- pnorm(Em1)
    prob2 <- pnorm(Em2)
    
    out <- cbind(m1 * prob1 + (1 - m1) * (1 - prob1),
          m2 * prob2 + (1 - m2) * (1 - prob2) )
  }
  out
}

pm1aw <- function(m1, a, w) {
  pz(1, a, w) * (
    apply(pm(m1, m2=1, 1, a, w), 1, prod) + 
      apply(pm(m1, m2=0, 1, a, w), 1, prod) )  +
    pz(0, a, w) * (
      apply(pm(m1, m2=1, 0, a, w), 1, prod) + 
        apply(pm(m1, m2=0, 0, a, w), 1, prod) ) 
}

pm2aw <- function(m2, a, w) {
  pz(1, a, w) * (
    apply(pm(m1=1, m2, 1, a, w), 1, prod) + 
      apply(pm(m1=0, m2, 1, a, w), 1, prod) )  +
    pz(0, a, w) * (
      apply(pm(m1=1, m2, 0, a, w), 1, prod) + 
        apply(pm(m1=0, m2, 0, a, w), 1, prod) ) 
}

pm1m2aw <- function(m1, m2, a1, a2, w) {
    pm1aw(m1, a1, w) * pm2aw(m2, a2, w)
}

pmaw <- function(m1, m2, a, w) {
  pz(1, a, w) * apply(pm(m1, m2, 1, a, w), 1, prod) +
    pz(0, a, w) * apply(pm(m1, m2, 0, a, w), 1, prod)
}


m_ratio <- function(m1, m2, z, a, w) {
  pm1aw(m1, a, w) * pm2aw(m2, a, w) / apply(pm(m1, m2, z, a, w), 1, prod)
}



my <- function(m1, m2, z, a, w, binary = TRUE) {
  mzaw <- 0.3*m1 + 0.3*m2 + 0.3*z + 0.3*a + rowSums(sqrt(0.13/ncol(w)) * w) # no interaction
  if (binary) {
    out <- pnorm(mzaw)
  }
  if (!binary) {
    out <- mzaw + rnorm(nrow(w), sd = sqrt(1-0.13))
  }
  out
}

u <- function(z, w, aprime, astar) {
    my(1, 1, z, aprime, w) * pmaw(1, 1, astar, w) +
        my(0, 1, z, aprime, w) * pmaw(0, 1, astar, w) +
        my(0, 0, z, aprime, w) * pmaw(0, 0, astar, w) +
        my(1, 0, z, aprime, w) * pmaw(1, 0, astar, w)
}

intu <- function(w, aprime, astar) {
    u(1, w, aprime, astar) * pz(1, aprime, w) +
      u(0, w, aprime, astar) * pz(0, aprime, w) 
}

intv <- function(m1, m2, w, aprime) {
    my(m1, m2, 1, aprime, w) * pz(1, aprime, w) +
        my(m1, m2, 0, aprime, w) * pz(0, aprime, w) 
}



gendata <- function(N = 1000) {
  # w1 <- rbinom(N, 1, .4)
  # w2 <- rbinom(N, 1, .4)
  w <- w_n(N, 3)
  
  a <- rbinom(N, 1, g(1, w))
  z <- rbinom(N, 1, pz(1, a, w))
  pm_jo <- pm(1, 1, z, a, w)
  m1 <- rbinom(N, 1, pm_jo[, 1])
  m2 <- rbinom(N, 1, pm_jo[, 2])
  y <- rbinom(N, 1, my(m1, m2, z, a, w))
  data.frame(W1 = w, A = a, Z1 = z, M1 = m1, M2 = m2, Y = y)
}




truth <- function(Nsuper = 1e4) {
  w <- w_n(Nsuper, 3)
  
  
  true_thetas <- rep(NA, 6)
  names(true_thetas) <- c("Y(0,M(t_jo=0))", "Y(1,M(t_jo=0))", "Y(1,M(t_jo=1))", 
                          "Y(1,M1(1),M2(1))", "Y(1,M1(0),M2(1))", "Y(1,M1(0),M2(0))")
  
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
      
      true_v <- intv(1, 1, w, t0) * pmaw(1, 1, t_jo, w) +
        intv(1, 0, w, t0) * pmaw(1, 0, t_jo, w) +
        intv(0, 1, w, t0) * pmaw(0, 1, t_jo, w) +
        intv(0, 0, w, t0) * pmaw(0, 0, t_jo, w)
      
      vname <- paste0("Y(", t0 ,",M(t_jo=", t_jo, "))")
      true_thetas[vname] <- mean(true_v)
    }
    if (is.na(param[["t_jo"]])) {
      t_1 <- param["t_1"]
      t_2 <- param["t_2"]
      
      true_v <- intv(1, 1, w, t0) * pm1m2aw(1, 1, t_1, t_2, w) +
        intv(1, 0, w, t0) * pm1m2aw(1, 0, t_1, t_2, w) +
        intv(0, 1, w, t0) * pm1m2aw(0, 1, t_1, t_2, w) +
        intv(0, 0, w, t0) * pm1m2aw(0, 0, t_1, t_2, w)
      
      vname <- paste0("Y(", t0 ,",M1(", t_1, "),", "M2(", t_2, ")", ")")
      # true_thetas[vname] <- weighted.mean(true_v, prob_w)
      true_thetas[vname] <- mean(true_v)
    }
  }
  
  true_effs <- effect(true_thetas)
  
  true_effs
}


trueval <- truth()

If <- function(dat) {
  w <- dat[, Wnames, drop = F]
  true_Ifs <- rep(NA, 6)
  names(true_Ifs) <- c("Y(0,M(t_jo=0))", "Y(1,M(t_jo=0))", "Y(1,M(t_jo=1))", 
                          "Y(1,M1(1),M2(1))", "Y(1,M1(0),M2(1))", "Y(1,M1(0),M2(0))")
  
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
      
      # EIF calculation
      ipwy <- (dat$A == t0) / g(t0, w)
      hmjo <- pmaw(dat$M1, dat$M2, t_jo, w) / pm(dat$M1, dat$M2, dat$Z1, t0, w)
      eify <- ipwy * hmjo / mean(ipwy * hmjo) * (dat$Y - my(dat$M1, dat$M2, dat$Z1, t0, w))
      
      ipwz <- ipwy
      eifz <- ipwz / mean(ipwz) * (u(dat$Z, w, t0, t_jo) - intu(w, t0, t_jo))
      
      ipwm <- (dat$A == t_jo) / g(t_jo, w)
      vbar <- intv(1, 1, w, t0) * pmaw(1, 1, t_jo, w) + intv(1, 0, w, t0) * pmaw(1, 0, t_jo, w) +
        intv(0, 1, w, t0) * pmaw(0, 1, t_jo, w) + intv(0, 0, w, t0) * pmaw(0, 0, t_jo, w)
      eifm <- ipwm / mean(ipwm) * (intv(dat$M1, dat$M2, w, t0) - vbar)
      
      # eif <- rescale_y(eify + eifz + eifm + mubarZ[, 1], bounds$bounds)
      eif <- eify + eifz + eifm + vbar
      theta <- mean(eif)
      
      vname <- paste0("Y(", t0 ,",M(t_jo=", t_jo, "))")
      vvbar[[vname]] <- mubarZ[, 1]
      
      thetas[vname] <- theta
      eifs[[vname]] <- eif
      
    }
    if (is.na(param[["t_jo"]])) {
      t_1 <- param["t_1"]
      t_2 <- param["t_2"]
      
      
      
      vname <- paste0("Y(", t0 ,",M1(", t_1, "),", "M2(", t_2, ")", ")")
      # true_thetas[vname] <- weighted.mean(true_v, prob_w)
      true_thetas[vname] <- mean(true_v)
    }
  }
  
}

