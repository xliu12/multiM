h_mjo <- function(h_z, g, em1m2, t_jo, t_obs) {
    `g(a'|w)` <- g[, gl("g({t_obs}|w)")]
    `g(a*|w)` <- g[, gl("g({t_jo}|w)")]
    `e(a'|m,w)` <- em1m2[, gl("e({t_obs}|m1,m2,w)")]
    `e(a*|m,w)` <- em1m2[, gl("e({t_jo}|m1,m2,w)")]
    h_z[, gl("h_z({t_obs})")] * `g(a'|w)` / `g(a*|w)` * `e(a*|m,w)` / `e(a'|m,w)`
}



h_m1 <- function(h_z1, g, em1, t_1, t_obs) {
    `g(a'|w)` <- g[, gl("g({t_obs}|w)")]
    `g(a*|w)` <- g[, gl("g({t_1}|w)")]
    `e(a'|m,w)` <- em1[, gl("e({t_obs}|m1,w)")]
    `e(a*|m,w)` <- em1[, gl("e({t_1}|m1,w)")]
    h_z1[, gl("h_z1({t_obs})")] * `g(a'|w)` / `g(a*|w)` * `e(a*|m,w)` / `e(a'|m,w)`
}


h_m2sec <- function(h_zm1, g, em2, t_2, t_obs) {
    `g(a'|w)` <- g[, gl("g({t_obs}|w)")]
    `g(a*|w)` <- g[, gl("g({t_2}|w)")]
    `e(a'|m,w)` <- em2[, gl("e({t_obs}|m2,w)")]
    `e(a*|m,w)` <- em2[, gl("e({t_2}|m2,w)")]
    h_zm1[, gl("h_zm1({t_obs})")] * `g(a'|w)` / `g(a*|w)` * `e(a*|m,w)` / `e(a'|m,w)`
}

h_m1sec <- function(h_zm2, g, em1, t_1, t_obs) {
    `g(a'|w)` <- g[, gl("g({t_obs}|w)")]
    `g(a*|w)` <- g[, gl("g({t_1}|w)")]
    `e(a'|m,w)` <- em1[, gl("e({t_obs}|m1,w)")]
    `e(a*|m,w)` <- em1[, gl("e({t_1}|m1,w)")]
    h_zm2[, gl("h_zm2({t_obs})")] * `g(a'|w)` / `g(a*|w)` * `e(a*|m,w)` / `e(a'|m,w)`
}
