#' matrix_to_estim
#'
#' @description
#' Function to turn a list of matrices into a vector of estimates/SEs. Performs
#' the inverse of estim_to_cov
#'
#' @param matrices_list List of diagonal Pk matrices and Sigmak matrices
#'
#' @return Vector of estimates/SEs
#'
#' @keywords internal
#'

matrix_to_estim <- function(matrices_list) {
  K <- length(matrices_list)/2
  P_vec <- unlist(lapply(matrices_list[1:K], diag))
  Sigma_vec <- unlist(lapply(matrices_list[(K+1):(2*K)], function(x){
    return(x[lower.tri(x, diag = TRUE)])
  }))
  return(c(P_vec, Sigma_vec))
}
