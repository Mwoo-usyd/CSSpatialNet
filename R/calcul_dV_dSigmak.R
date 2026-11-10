calcul_dV_dSigmak <- function(Ainv_list,
                              D_list,
                              n) {
  K <- length(Ainv_list)
  p <- nrow(Ainv_list[[1]])/n
  dV_dPk_list <- lapply(c(1:K), function(k) {
    ainv_d <- t(t(Ainv_list[[k]]) * diag(D_list[[k]]))
    d_ainv <- Ainv_list[[k]] * diag(D_list[[k]])
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
