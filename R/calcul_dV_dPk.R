#' calcul_dV_dPk
#'
#' @description
#' Function to calculate the derivative of V (covariance matrix of vec(X)) with
#' respect to each spatial autocorrelation parameter
#'
#' @param Ainv_list List of matrices of the form
#' \eqn{(I_{n} - P_{k,i} \otimes W)^{-1}} where \eqn{P_{k,i}} is the spatial
#' autocorrelation of the i'th gene in the k'th cell type
#' @param D_list List of expanded proportions of length K. Each is a diagonal
#' \eqn{np \times np} matrix whose values are the k'th cell type proportions,
#' repeated p times.
#' @param Sigma_list List of cell-type-specific covariance matrices
#' @param W Spatial weight matrix
#'
#' @return List of matrices, each corresponding to the derivative of V with
#' respect to a spatial autocorrelation parameter
#'
#' @importFrom Matrix Matrix
#' @importFrom purrr list_flatten
#'
#' @keywords internal
#'

calcul_dV_dPk <- function(Ainv_list,
                          D_list,
                          Sigma_list,
                          W) {
  K <- length(Ainv_list)
  n <- nrow(W)
  p <- nrow(Ainv_list[[1]])/n
  dV_dPk_list <- lapply(c(1:K), function(k) {
    ainv_d <- t(t(as.matrix(Ainv_list[[k]])) * diag(D_list[[k]]))
    d_ainv <- Ainv_list[[k]] * diag(D_list[[k]])

    sigmak_kron_In <- Matrix(kronecker(Sigma_list[[k]],
                                       diag(n)),
                             sparse = TRUE)

    grad_list <- list()
    for (i in 1:p) {
      mat1 <- matrix(rep(0, p^2), nrow = p)
      mat1[i,i] <- 1

      kron_m1_W <- Matrix(kronecker(mat1, W), sparse = TRUE)

      inside <- kron_m1_W %*% Ainv_list[[k]] %*% sigmak_kron_In +
        sigmak_kron_In %*% Ainv_list[[k]] %*% kron_m1_W

      deriv_current <- d_ainv %*% inside %*% ainv_d
      grad_list <- grad_list %>% append(list(deriv_current))
    }
    return(grad_list)
  })
  return(list_flatten(dV_dPk_list))
}
