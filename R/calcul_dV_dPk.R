calcul_dV_dPk <- function(Ainv_list, D_list, Sigma_list, W) {
  K <- length(Ainv_list)
  n <- nrow(W)
  p <- nrow(Ainv_list[[1]])/n
  dV_dPk_list <- lapply(c(1:K), function(k) {
    ainv_d <- t(t(Ainv_list[[k]]) * diag(D_list[[k]]))
    d_ainv <- Ainv_list[[k]] * diag(D_list[[k]])

    sigmak_kron_In <- Matrix(kronecker(Sigma_list[[k]], diag(n)), sparse = TRUE)

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
