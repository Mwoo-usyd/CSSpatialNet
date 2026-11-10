calcul_d2V_dPk_dPk <- function(Ainv_list, D_Ainv_list, Ainv_D_list, W, Sigma_list) {
  K <- length(Ainv_list)
  n <- nrow(W)
  p <- nrow(Sigma_list[[1]])

  grads_over_k <- lapply(c(1:K), function(k) {
    sigmak_kron_In <- Matrix(kronecker(Sigma_list[[k]], diag(n)), sparse = TRUE)
    grads_over_dim1 <- lapply(c(1:p), function(i){
      mat1i <- Matrix(0, nrow = p, ncol = p, sparse = TRUE)
      mat1i[i,i] <- 1
      i_component <- Matrix(kronecker(mat1i, W), sparse = TRUE)
      sigmak_i_component <- sigmak_kron_In %*% Ainv_list[[k]] %*% i_component
      i_sigmak_component <- i_component %*% Ainv_list[[k]] %*% sigmak_kron_In
      grads_over_dim2 <- lapply(c(1:i), function(j){
        mat1j <- Matrix(0, nrow = p, ncol = p, sparse = TRUE)
        mat1j[j,j] <- 1
        j_component <- Matrix(kronecker(mat1j, W), sparse = TRUE)
        j_ainv <- j_component %*% Ainv_list[[k]]
        ainv_j <- Ainv_list[[k]] %*% j_component

        # create 6 combinations
        i_sigmak_j <- i_sigmak_component %*% ainv_j
        sigmak_i_j <- sigmak_i_component %*% ainv_j
        j_i_sigmak <- j_ainv %*% i_sigmak_component
        j_sigmak_i <- j_ainv %*% sigmak_i_component
        i_j_sigmak <- i_component %*% Ainv_list[[k]] %*% j_ainv %*% sigmak_kron_In
        sigmak_j_i <- sigmak_kron_In %*% Ainv_list[[k]] %*% j_ainv %*% i_component
        sum6 <- i_sigmak_j + sigmak_i_j + j_i_sigmak + j_sigmak_i + i_j_sigmak + sigmak_j_i
        final_grad <- D_Ainv_list[[k]] %*% sum6 %*% Ainv_D_list[[k]]
        return(final_grad)
      })
      return(list_flatten(grads_over_dim2))
    })
    return(list_flatten(grads_over_dim1))
  })
  return(list_flatten(grads_over_k))
}
