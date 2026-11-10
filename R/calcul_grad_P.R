# function to calculate the gradient of log-likelihood with respect to all parameters (Pk diagonal and Sigmak diagonal and off-diagonal)
calcul_grad_P <- function(current_P_list, # list of K P matrices, each pxp
                          current_Sigma_list, # list of K Sigma matrices, each pxp
                          D_list,
                          W,
                          svd_obj, # created from svd(W), outside of this so as to not redo it many times. Includes slots v_inv and u_inv
                          x, # n * p matrix
                          Ainv_list = NULL,
                          D_Ainv_list = NULL,
                          Ainv_D_list = NULL,
                          current_V = NULL) {
  K <- length(current_Sigma_list)
  n <- length(svd_obj$d)
  p <- nrow(current_P_list[[1]])

  if (is.null(Ainv_list)) {
    Ainv_list <- get_ainv(current_P_list,
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
  if (is.null(current_V)) {
    current_V <- calcul_V_P(current_Ainv_list = Ainv_list,
                            current_Sigma_list = current_Sigma_list,
                            D_list = D_list)
  }

  term1 <- (1/2) * current_V - 0.5 * tcrossprod(as.matrix(c(x))) # this is not dependent on i or if we are looking at rho or sigma2 # here want it to be np*np
  # V_inv <- chol2inv(chol(current_V)) # Up to twice as fast as solve(current_V)
  # V_inv <- ginv(as.matrix(current_V))
  # F_mat <- V_inv %*% term1 %*% V_inv # meant to be t(term1) but it is symmetric anyway
  intermediate_product <- solve(current_V, term1) # this is a faster way to get current_V^-1 %*% term1 %*% current_V^-1
  F_mat <- t(solve(t(current_V), t(intermediate_product)))

  grad_list <- list() # list of gradient matrices, length 2K since has parameters P_i and Sigma_i. Each element is a pxp matrix of derivatives

  G_list <- lapply(c(1:K), function(k){
    G_mat <- Ainv_D_list[[k]] %*% F_mat %*% D_Ainv_list[[k]]
    return(G_mat)
  })

  for (k in 1:K) { # first fill out first K terms with P_k derivatives
    current_P_i_deriv <- matrix(0, nrow = p, ncol = p)

    kron_component <- Matrix(kronecker(current_Sigma_list[[k]], diag(n)), sparse = TRUE)
    first_part <- Ainv_list[[k]] %*% kron_component
    second_part <- kron_component %*% Ainv_list[[k]]

    for (i in 1:p) {
      mat1 <- matrix(rep(0, p^2), nrow = p)
      mat1[i,i] <- 1

      kron_m1_W <- Matrix(kronecker(mat1, W), sparse = TRUE)

      inside <- kron_m1_W %*% first_part + second_part %*% kron_m1_W

      deriv_current <- -sum(diag(G_list[[k]] %*% inside))
      current_P_i_deriv[i,i] <- deriv_current
    }
    grad_list <- grad_list %>% append(list(current_P_i_deriv))
  }
  for (k in 1:K) { # next fill out remaining K terms with Sigma_i derivatives
    G_mat_current <- G_list[[k]]

    current_Sigma_i_deriv <- matrix(0, nrow = p, ncol = p)
    for (i in 1:p) {
      for (j in 1:p) {
        rows_to_take <- ((i-1)*n+1):(i*n)
        cols_to_take <- ((j-1)*n+1):(j*n)
        G_subset <- G_mat_current[rows_to_take, cols_to_take]
        current_Sigma_i_deriv[i,j] <- -sum(diag(G_subset))
      }
    }
    grad_list <- grad_list %>% append(list(current_Sigma_i_deriv))
  }
  return(grad_list)
}
