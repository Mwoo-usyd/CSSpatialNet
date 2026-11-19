#' calcul_dl_dVinv
#'
#' @description
#' Function to calculate the derivative of the log-likelihood with respect to
#' the inverse of the covariance matrix of vec(X)
#'
#' @param vec_X vectorised gene expression matrix
#' @param V (possible) covariance matrix of vec(X)
#'
#' @return Matrix of derivative of log-likelihood with respect to V-inverse
#'
#' @keywords internal

calcul_dl_dVinv <- function(vec_X,
                            V) {
  return((1/2) * V - 0.5 * tcrossprod(as.matrix(c(vec_X))))
}
