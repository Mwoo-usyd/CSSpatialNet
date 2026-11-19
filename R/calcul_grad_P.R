#' calcul_grad_P
#'
#' @description
#' Function to calculate the gradient of the log-likelihood with respect to
#' each spatial autocorrelation and cell-type-specific covariance matrix
#' parameter
#'
#' @param P_list List of cell-type specific diagonal matrices of spatial
#' autocorrelation
#' @param Sigma_list List of cell-type-specific covariance matrices
#' @param D_list List of expanded proportions of length K. Each is a diagonal
#' \eqn{np \times np} matrix whose values are the k'th cell type proportions,
#' repeated p times.
#' @param W Spatial weight matrix
#' @param svd_obj Output of performing SVD on the spatial weight matrix. Also
#' includes slots for the inverses of U and V, as well as t(U) * V.
#' @param vec_X Vectorised gene expression matrix
#' @param Ainv_list List of matrices of the form
#' \eqn{(I_{n} - P_{k,i} \otimes W)^{-1}} where \eqn{P_{k,i}} is the spatial
#' autocorrelation of the i'th gene in the k'th cell type
#' @param D_Ainv_list List (length K) of each
#' \eqn{D_k \times (I_{np}-P_{k}W)^{-1}}
#' @param Ainv_D_list List (length K) of each
#' \eqn{(I_{np}-P_{k}W)^{-1} \times D_k}
#' @param V Covariance matrix of vec(X)
#'
#' @return List of gradients of the log-likelihood with respect to the
#' (diagonal) matrices of cell-type-specific spatial autocorrelation, and the
#' gene-gene covariance matrices
#'
#' @importFrom Matrix Matrix
#'
#' @keywords internal

calcul_grad_P <- function(P_list,
                          Sigma_list,
                          D_list,
                          W,
                          svd_obj,
                          vec_X,
                          Ainv_list = NULL,
                          D_Ainv_list = NULL,
                          Ainv_D_list = NULL,
                          V = NULL) {
  K <- length(Sigma_list)
  n <- length(svd_obj$d)
  p <- nrow(P_list[[1]])

  if (is.null(Ainv_list)) {
    Ainv_list <- get_ainv(P_list,
                          d = svd_obj$d,
                          U_inv = svd_obj$u_inv,
                          V_inv = svd_obj$v_inv,
                          tUxV = svd_obj$tUxV)
  }
  if (is.null(Ainv_D_list)) {
    Ainv_D_list <- lapply(c(1:K), function(k){
      return(t(t(Ainv_list[[k]]) * diag(D_list[[k]])))
    })
  }
  if (is.null(D_Ainv_list)) {
    D_Ainv_list <- lapply(c(1:K), function(k){
      return(Ainv_list[[k]] * diag(D_list[[k]]))
    })
  }
  if (is.null(V)) {
    V <- calcul_V_P(Ainv_list = Ainv_list,
                    Sigma_list = Sigma_list,
                    D_list = D_list)
  }

  term1 <- (1/2) * V - 0.5 * tcrossprod(as.matrix(c(vec_X)))
  intermediate_product <- solve(V, as.matrix(term1))
  F_mat <- t(solve(t(as.matrix(V)), t(intermediate_product)))

  grad_list <- list()

  G_list <- lapply(c(1:K), function(k){
    G_mat <- Ainv_D_list[[k]] %*% F_mat %*% D_Ainv_list[[k]]
    return(G_mat)
  })

  for (k in 1:K) {
    current_P_i_deriv <- matrix(0, nrow = p, ncol = p)

    kron_component <- Matrix(kronecker(Sigma_list[[k]], diag(n)), sparse = TRUE)
    first_part <- Ainv_list[[k]] %*% kron_component
    second_part <- kron_component %*% Ainv_list[[k]]

    for (i in 1:p) {
      mat1 <- matrix(rep(0, p^2), nrow = p)
      mat1[i,i] <- 1

      kron_m1_W <- Matrix(kronecker(mat1, W), sparse = TRUE)

      inside <- kron_m1_W %*% first_part + second_part %*% kron_m1_W

      deriv_current <- -sum(diag(as.matrix(G_list[[k]] %*% inside)))
      current_P_i_deriv[i,i] <- deriv_current
    }
    grad_list <- grad_list %>% append(list(current_P_i_deriv))
  }
  for (k in 1:K) {
    G_mat_current <- G_list[[k]]

    current_Sigma_i_deriv <- matrix(0, nrow = p, ncol = p)
    for (i in 1:p) {
      for (j in 1:p) {
        rows_to_take <- ((i-1)*n+1):(i*n)
        cols_to_take <- ((j-1)*n+1):(j*n)
        G_subset <- G_mat_current[rows_to_take, cols_to_take]
        current_Sigma_i_deriv[i,j] <- -sum(diag(as.matrix(G_subset)))
      }
    }
    grad_list <- grad_list %>% append(list(current_Sigma_i_deriv))
  }
  return(grad_list)
}
