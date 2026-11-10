calcul_dl_dVinv <- function(z, current_V) {
  return((1/2) * current_V - 0.5 * tcrossprod(as.matrix(c(z))))
}
