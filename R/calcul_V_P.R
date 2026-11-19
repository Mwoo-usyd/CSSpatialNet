
#' calcul_V_P
#'
#' @description
#' Function to calculate the covariance of vec(X), which is the vectorised
#' spot-level gene expression data, given some parameters and proportions.
#'
#' @param Ainv_list List of matrices of the form
#' \eqn{(I_{n} - P_{k,i} \otimes W)^{-1}} where \eqn{P_{k,i}} is the spatial
#' autocorrelation of the i'th gene in the k'th cell type
#' @param Sigma_list List of cell-type-specific covariance matrices
#' @param D_list List of expanded proportions of length K. Each is a diagonal
#' \eqn{np \times np} matrix whose values are the k'th cell type proportions,
#' repeated p times.
#'
#' @return Covariance matrix of vectorised spot-level gene expression matrix
#' using the input parameters. Dimension \eqn{np \times np}
#'
#' @keywords internal

calcul_V_P <- function(Ainv_list,
                       Sigma_list,
                       D_list) {
  K <- length(Sigma_list)
  p <- nrow(Sigma_list[[1]])
  n <- nrow(D_list[[1]])/p
  V <- matrix(0, nrow = n*p, ncol = n*p)
  for (i in 1:K) {
    middle <- Ainv_list[[i]] %*%
      kronecker(Sigma_list[[i]], diag(n)) %*%
      Ainv_list[[i]]
    current_add <- diag(D_list[[i]]) *
      middle *
      rep(diag(D_list[[i]]), each = n*p)
    V <- V + current_add
  }
  return(V)
}
