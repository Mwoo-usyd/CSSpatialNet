#' calcul_dV_dSigmak
#'
#' @description
#' Function to calculate the derivative of V (covariance matrix of vec(X)) with
#' respect to each parameter in the covariance matrix of each cell type's gene
#' expression matrix.
#'
#' @param Ainv_list List of matrices of the form
#' \eqn{(I_{n} - P_{k,i} \otimes W)^{-1}} where \eqn{P_{k,i}} is the spatial
#' autocorrelation of the i'th gene in the k'th cell type
#' @param D_list List of expanded proportions of length K. Each is a diagonal
#' \eqn{np \times np} matrix whose values are the k'th cell type proportions,
#' repeated p times.
#' @param n Sample size
#'
#' @return List of matrices, each corresponding to the derivative of V with
#' respect to a gene-gene covariance matrix parameter
#'
#' @importFrom Matrix Matrix
#' @importFrom purrr list_flatten
#'
#' @keywords internal
#'

calcul_dV_dSigmak <- function(Ainv_list,
                              D_list,
                              n) {
  K <- length(Ainv_list)
  p <- nrow(Ainv_list[[1]])/n
  dV_dPk_list <- lapply(c(1:K), function(k) {
    ainv_d <- t(t(as.matrix(Ainv_list[[k]])) * diag(as.matrix(D_list[[k]])))
    d_ainv <- Ainv_list[[k]] * diag(as.matrix(D_list[[k]]))
    grad_list <- list()
    for (i in 1:p) {
      for (j in i:p) {
        mat1 <- matrix(rep(0, p^2), nrow = p)
        mat1[j,i] <- 1
        mat1[i,j] <- 1

        kron_m1_In <- Matrix(kronecker(mat1, diag(n)), sparse = TRUE)
        deriv_current <- d_ainv %*% kron_m1_In %*% ainv_d
        grad_list <- grad_list %>% append(list(deriv_current))
      }
    }
    return(grad_list)
  })
  return(list_flatten(dV_dPk_list))
}
