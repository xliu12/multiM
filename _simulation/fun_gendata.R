pa <- function(a, given, gen_a) {
  prob1 <- pnorm(given, mean = 0, sd = sqrt(1 - gen_a$icca))
  a * prob1 + (1 - a) * (1 - prob1)
}

pm1 <- function(m1, a, given, gen_m) {
  latent <- gen_m[["m1_on_a"]] * a + given
  prob1 <- pnorm(latent, mean = 0, sd = sqrt(1 - gen_m[["iccm1"]]))
  m1 * prob1 + (1 - m1) * (1 - prob1)
}

pm2 <- function(m2, m1, a, given, gen_m) {
  latent <- gen_m[["m2_on_m1"]] * m1 + gen_m[["m2_on_a"]] * a +
    gen_m[["m2_on_am1"]] * a * m1 + 
    given
  
  prob1 <- pnorm(latent, mean = 0, sd = sqrt(1 - gen_m[["iccm2"]]))
  m2 * prob1 + (1 - m2) * (1 - prob1)
}

pm2a <- function(m2, a, given, gen_m) {
    pm2(m2, 1, a, given, gen_m) * pm1(1, a, given, gen_m) +
        pm2(m2, 0, a, given, gen_m) * pm1(0, a, given, gen_m)
}


my <- function(m2, m1, a, given, gen_y, binary = TRUE) {
  latent <- gen_y[["y_on_m2"]] * m2 + gen_y[["y_on_m1"]] * m1 + gen_y[["y_on_a"]] * a + 
    gen_y[["y_on_am2"]] * a * m2 + gen_y[["y_on_am1"]] * a * m1 + gen_y[["y_on_m1m2"]] * m1 * m2 + 
    gen_y[["y_on_am1m2"]] * a * m1 * m2 + 
    given
  
  if (binary) {
    cond_mean <- pnorm(latent, mean = 0, sd = sqrt(1 - gen_y[["iccy"]]))
  }
  if (!binary) {
    cond_mean <- plogis(latent) # i.e., expit(.) = exp(.)/(1+exp(.))
  }
  
  cond_mean
}

# u <- function(z, w, aprime, astar) {
#     my(1, z, aprime, w) * pmaw(1, astar, w) + my(0, z, aprime, w) * pmaw(0, astar, w)
# }
# 
# intu <- function(w, aprime, astar) {
#     u(1, w, aprime, astar) * pz(1, aprime, w) +
#         u(0, w, aprime, astar) * pz(0, aprime, w)
# }
# 
# intv <- function(m, w, aprime) {
#     my(m, 1, aprime, w) * pz(1, aprime, w) +
#         my(m, 0, aprime, w) * pz(0, aprime, w)
# }
# 
# gendata <- function(N) {
#     w1 <- rbinom(N, 1, .4)
#     w <- data.frame(W1 = w1)
#     a <- rbinom(N, 1, g(1))
#     z <- rbinom(N, 1, pz(1, a, w))
#     m <- rbinom(N, 1, pm(1, z, a, w))
#     y <- rbinom(N, 1, my(m, z, a, w))
#     data.frame(W1 = w1, A = a, Z = z, M = m, Y = y)
# }
# 
# # h <- pmaw(m, astar, w) / pm(m, z, aprime, w)
# 
# If <- function(dat, aprime, astar) {
#     w <- dat[, "W1", drop = F]
#     
#     ipwy <- (dat$A == aprime) / g(aprime)
#     hm <- pmaw(dat$M, astar, w) / pm(dat$M, dat$Z, aprime, w)
#     eify <- ipwy * hm * (dat$Y - my(dat$M, dat$Z, aprime, w))
#     
#     ipwz <- ipwy
#     eifz <- ipwz * (u(dat$Z, w, aprime, astar) - intu(w, aprime, astar))
#     
#     ipwm <- (dat$A == astar) / g(astar)
#     vbar <- intv(1, w, aprime) * pmaw(1, astar, w) + intv(0, w, aprime) * pmaw(0, astar, w)
#     eifm <- ipwm * (intv(dat$M, w, aprime) - vbar)
#     
#     eify + eifz + eifm + vbar
# }
# 
# truth <- function() {
#     w <- expand.grid(W1 = c(0, 1))
# 
#     prob_w <- vector("numeric", nrow(w))
#     for (i in 1:nrow(w)) {
#         w1 <- w[i, "W1"] * 0.4 + (1 - w[i, "W1"]) * 0.6
#         prob_w[i] <- w1
#     }
# 
#     aprime <- astar <- 1
#     v_11 <- intv(1, w, aprime) * pmaw(1, astar, w) + intv(0, w, aprime) * pmaw(0, astar, w)
# 
#     astar <- 0
#     v_10 <- intv(1, w, aprime) * pmaw(1, astar, w) + intv(0, w, aprime) * pmaw(0, astar, w)
# 
#     aprime <- 0
#     v_00 <- intv(1, w, aprime) * pmaw(1, astar, w) + intv(0, w, aprime) * pmaw(0, astar, w)
# 
#     c("11" = weighted.mean(v_11, prob_w),
#       "10" = weighted.mean(v_10, prob_w),
#       "00" = weighted.mean(v_00, prob_w),
#       "indirect" = weighted.mean(v_11 - v_10, prob_w),
#       "direct" = weighted.mean(v_10 - v_00, prob_w))
# }
