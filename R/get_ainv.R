# function to speed up calculation of (I - Pk kron W)^{-1}, using the fact tha Pk is diagonal
get_ainv <- function(P_list, U_inv, d, V_inv, tUxV) { # where U, d, V are from svd(W), and _inv is inverse. d is a vector
  n <- length(d)
  p <- nrow(P_list[[1]])

  Ainv_all_K <- lapply(P_list, function(Pk) {
    diag_vals <- diag(Pk) # extract diagonal values (off-diags are all 0 anyway)
    each_gene_inv <- lapply(diag_vals, function(Pki) { # apply function to each scalar Pk[i,i]
      # result <- t(V_inv) %*% solve(tUxV - diag(Pki * d)) %*% U_inv
      result <- crossprod(V_inv, solve(tUxV - diag(Pki * d))) %*% U_inv
      return(result)
    })
    return(bdiag(each_gene_inv))
  })
  return(Ainv_all_K)
}
