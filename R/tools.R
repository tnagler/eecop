compute_pseudo_obs <- function(x, margins) {
  if (is.numeric(x) && all(dim(x) == c(length(margins), 1)))
    x <- t(x)
  d <- NCOL(x)
  u <- lapply(seq_len(d), function(k) margins[[k]](x[, k]))
  u <- do.call(cbind, u)

  i_fct <- which(sapply(x, is.factor))
  if (length(i_fct) > 0) {
    for (k in i_fct) {
      lv <- as.numeric(x[, k]) - 1
      lv0 <- which(lv == 0)
      lv[lv0] <- 1
      xlv <- factor(levels(x[, k])[lv], levels = levels(x[, k]))
      u_sub <- margins[[k]](xlv)
      u_sub[lv0] <- 0
      u <- cbind(u, u_sub)
    }
  }
  cut_01(u)
}


cut_01 <- function(x, gap = 1e-10) {
  pmin(pmax(x, gap), 1 - gap)
}

combine_margins <- function(V, U, q, p) {
  cbind(V[, seq_len(q)], U[, seq_len(p)], V[, -seq_len(q)], U[, -seq_len(p)])
}

