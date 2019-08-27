#' @importFrom stats qnorm
#' @noRd
fit_copula_gaussian <- function(u, weights) {
    rho <- wdm::wdm(qnorm(u), weights = weights)
    cop_model <- copula::normalCopula(param = copula::P2p(rho),
                                      dim = ncol(u),
                                      dispstr = "un")
    function(u_new) copula::dCopula(u_new, cop_model)
}

#' @noRd
fit_copula_vine <- function(u, weights, ...) {
    cop_model <- rvinecopulib::vinecop(u, weights = weights, ...)
    function(u_new) rvinecopulib::dvinecop(u_new, cop_model)
}

#' @importFrom stats qnorm dnorm
#' @noRd
fit_copula_kde <- function(u, weights) {
    x <- qnorm(u)
    bw <- nrow(x)^(-2 / (ncol(x) + 2)) * cov(x)
    function(u_new) {
        x_new <- qnorm(u_new)
        n <- nrow(x)
        if (length(weights) == 0) {
            weights <- rep(1, n)
        }
        bws <- replicate(n, bw, simplify = FALSE)
        bws <- do.call(rbind, bws)
        f <- ks::dmvnorm.mixt(
            x = x_new,
            mus = x,
            Sigmas = bws,
            props = weights / sum(weights)
        )
        f / apply(dnorm(x_norm), 1, prod)
    }
}

#' @noRd
fit_copula_bernstein <- function(u, weights) {
    m <- ceiling(mult * nrow(u)^(2 / (4 + ncol(u))))
    fit <- fit_bern_coefs(u, m)
    function(u_new) eval_bern(u_new, fit)
}

fit_bern_coefs <- function(u, m, weights = numeric(0)) {
    if (!length(weights))
        weights <- rep(1, NROW(u))
    d <- NCOL(u)
    ucut <- lapply(seq_len(d), function(j) cut(u[, j], seq(0, 1, len = m + 1)))
    joined <- do.call(paste, ucut)
    lvls <- do.call(expand.grid, lapply(ucut, levels))
    lvls <- apply(lvls, 1, paste, collapse = " ")
    joined <- factor(joined, levels = lvls)
    coef <- tapply(weights, joined, sum) / NROW(u)
    basis_ind <- expand.grid(replicate(d, seq_len(m + 1), simplify = FALSE)) - 1
    ind <- which(cf_arr > 0)
    list(basis_ind = as.matrix(basis_ind[ind, ]),
         coef = coef[ind],
         m = m,
         d = d)
}

eval_bern <- function(u, fit) {
    c <- 0
    for (k in seq_len(nrow(fit$basis_ind))) {
        p <- fit$coef[k]
        for (j in seq_len(d)) {
            p <- p * bern_poly(u[, j], fit$basis_ind[k, j], fit$m)
        }
        c <- c + p
    }
    c
}

