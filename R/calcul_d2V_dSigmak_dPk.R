#' calcul_d2V_dSigmak_dPk
#'
#' @description
#' Calculates second derivative of V (covariance matrix of vec(X)) with respect
#' to a spatial autocorrelation parameter and a gene variance or covariance
#' parameter (both belonging to the same cell type)
#'
#' @param Ainv_list List of matrices of the form
#' \eqn{(I_{n} - P_{k,i} \otimes W)^{-1}} where \eqn{P_{k,i}} is the spatial
#' autocorrelation of the i'th gene in the k'th cell type
#' @param D_list List of expanded proportions of length K. Each is a diagonal
#' \eqn{np \times np} matrix whose values are the k'th cell type proportions,
#' repeated p times.
#' @param W Spatial weight matrix
#' @param D_Ainv_list List (length K) of each
#' \eqn{D_k \times (I_{np}-P_{k}W)^{-1}}
#' @param Ainv_D_list List (length K) of each
#' \eqn{(I_{np}-P_{k}W)^{-1} \times D_k}
#'
#' @return List of second derivative matrices, each corresponding to a
#' unique combination of cell type, gene, and covariance matrix entry, where the
#' matrix is the second derivative of V with respect to the gene and the
#' covariance matrix entry within one cell type.
#'
#' @importFrom Matrix Matrix
#' @importFrom purrr list_flatten
#'
#' @keywords internal

calcul_d2V_dSigmak_dPk <- function(Ainv_list,
                                   D_list,
                                   W,
                                   Ainv_D_list,
                                   D_Ainv_list) {
  K <- length(Ainv_list)
  n <- nrow(W)
  p <- nrow(Ainv_list[[1]])/n

  grads_over_K <- lapply(c(1:K), function(k){
    ainv_d <- Ainv_D_list[[k]]
    d_ainv <- D_Ainv_list[[k]]
    grads_over_Pk <- lapply(c(1:p), function(j){
      matj <- Matrix(0, nrow = p, ncol = p, sparse = TRUE)
      matj[j,j] <- 1
      W_component <- Matrix(kronecker(matj,W), sparse = TRUE)
      Dk_Ainvk_Wterm_Ainvk <- d_ainv %*% W_component %*% Ainv_list[[k]]
      Ainvk_Wterm_Ainvk_Dk <- Ainv_list[[k]] %*% W_component %*% ainv_d
      grads_over_dim1 <- lapply(c(1:p), function(i1){
        grads_over_dim2 <- lapply(c(1:i1), function(i2){
          mat1 <- Matrix(0, nrow = p, ncol = p, sparse = TRUE)
          mat1[i1,i2] <- 1
          mat1[i2,i1] <- 1
          In_component <- kronecker(mat1, diag(n))
          term1 <- Dk_Ainvk_Wterm_Ainvk %*% In_component %*% ainv_d
          term2 <- d_ainv %*% In_component %*% Ainvk_Wterm_Ainvk_Dk
          return(term1 + term2)
        })
        return(list_flatten(grads_over_dim2))
      })
      return(list_flatten(grads_over_dim1))
    })
    return(list_flatten(grads_over_Pk))
  })
  return(list_flatten(grads_over_K))
}
