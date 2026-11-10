# function to calculate the covariance matrix of vec(Z) given Sigmak, Ainvk (calculated from Pk, W), and the proportions list. Has P in the name to represent that this is for p > 1
calcul_V_P <- function(Ainv_list,
                       Sigma_list,
                       D_list) { # calculate the covariance matrix for currently-attempted values of all parameters
  K <- length(Sigma_list)
  p <- nrow(Sigma_list[[1]])
  n <- nrow(D_list[[1]])/p
  V <- matrix(0, nrow = n*p, ncol = n*p)
  for (i in 1:K) {
    middle <- Ainv_list[[i]] %*%
      kronecker(Sigma_list[[i]], diag(n)) %*%
      Ainv_list[[i]]
    current_add <- diag(D_list[[i]]) * middle * rep(diag(D_list[[i]]), each = n*p)
    V <- V + current_add
  }
  return(V)
}
