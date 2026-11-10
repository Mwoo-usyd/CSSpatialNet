# function to calculate the covariance matrix of vec(Z) given Sigmak, Ainvk (calculated from Pk, W), and the proportions list. Has P in the name to represent that this is for p > 1
calcul_V_P <- function(current_Ainv_list,
                       current_Sigma_list,
                       D_list) { # calculate the covariance matrix for currently-attempted values of all parameters
  K <- length(current_Sigma_list)
  p <- nrow(current_Sigma_list[[1]])
  n <- nrow(D_list[[1]])/p
  current_V <- matrix(0, nrow = n*p, ncol = n*p)
  for (i in 1:K) {
    middle <- current_Ainv_list[[i]] %*%
      kronecker(current_Sigma_list[[i]], diag(n)) %*%
      current_Ainv_list[[i]]
    current_add <- diag(D_list[[i]]) * middle * rep(diag(D_list[[i]]), each = n*p)
    current_V <- current_V + current_add
  }
  return(current_V)
}
