# Tianchen Qian
# 2021.05.30

# The following "estimate" function is copied from xgeepack.R,
# and I added the degrees of freedom in the output table.
# Please do not change it.
estimate <- function(x, combos = NULL, omnibus = FALSE, null = 0,
                     small = TRUE, conf.int = 0.95, normal = FALSE, ...) {
    if (is.null(combos)) {
        combos <- diag(length(coef(x)))
        rownames(combos) <- names(coef(x))
        omnibus <- FALSE
    }
    est <- combos %*% coef(x)
    if (nrow(est) != length(null)) null <- rep(null[1], nrow(est))
    ## apply Mancl and DeRouen's (2001) small sample correction
    if (is.logical(small)) small <- small * 50
    n <- cluster.number(x, overall = FALSE)
    d1 <- if (omnibus) nrow(combos)
    else apply(combos != 0, 1, sum)
    d2 <- n - length(coef(x))
    ## apply Hotelling's T-squared test, following Liao et al. (2016)
    if (n <= small & !normal) {
        type <- "Hotelling"
        adj <- d1 * (d1 + d2 - 1) / d2
        qfun <- function(p) mapply(qf, p = p, df1 = d1, df2 = d2) / adj
        pfun <- function(q) 1 - mapply(pf, q = q * adj, df1 = d1, df2 = d2)
    }
    else {
        type <- "Wald"
        qfun <- if (normal) function(p) qnorm((1 + p) / 2)
        else function(p) mapply(qf, p = p, df1 = d1, df2 = d2)
        pfun <- if (normal) function(q) 1 - mapply(pchisq, q = q, df = d1)
        else function(q) 1 - mapply(pf, q = q, df1 = d1, df2 = d2)
    }
    var.est <- combos %*% vcov(x, small = small, ...) %*% t(combos)
    se.est <- sqrt(diag(var.est))
    crit <- sqrt(qfun(conf.int))
    lcl <- est - se.est * crit
    ucl <- est + se.est * crit
    stat <- if (omnibus) rep(t(est - null) %*% solve(var.est) %*% (est - null), d1)
    else (est - null)^2 / diag(var.est)
    pvalue <- pfun(stat)
    out <- cbind(est, lcl, ucl, se.est, stat, d1, d2, pvalue)
    rownames(out) <- rownames(combos)
    colnames(out) <- c("Estimate",
                       paste0(round(conf.int * 100), "% ", c("LCL", "UCL")),
                       "SE", type, "df1", "df2", "p-value")
    class(out) <- c("estimate", "matrix")
    out
}