calcul_dl_dVinv <- function(vec_X,
                            V) {
  return((1/2) * V - 0.5 * tcrossprod(as.matrix(c(vec_X))))
}
