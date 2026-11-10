calcul_d2V_dSigmak_dPk <- function(Ainv_list, D_list, W, Ainv_D_list, D_Ainv_list) {
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
