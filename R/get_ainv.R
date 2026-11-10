# function to speed up calculation of (I - Pk kron W)^{-1}, using the fact tha Pk is diagonal

#' get_ainv
#'
#' @description
#' Function to calculate K matrices of the form (I_{np} - P_{k}*W)-1
#'
#' @param P_list List of cell-type specific diagonal matrices of
#' spatial autocorrelation
#' @param U_inv Inverse of U, which is an output of SVD(weight matrix)
#' @param d An output of SVD(weight matrix)
#' @param V_inv Inverse of V, which is an output of SVD(weight matrix)
#' @param tUxV transpose of U times V, where U and V are outputs
#' of SVD(weight matrix)
#'
#' @return List of K matrices of the form (I_{np} - P_{k}*W)-1
#'
#' @keywords internal

get_ainv <- function(P_list,
                     U_inv,
                     d,
                     V_inv,
                     tUxV) { # where U, d, V are from svd(W), and _inv is inverse. d is a vector
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
