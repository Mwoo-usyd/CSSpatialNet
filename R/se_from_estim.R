#' se_from_estim
#'
#' @description
#' Function to calculate the standard errors of each parameter through the
#' Fisher Information matrix.
#'
#' @param estim Matrix of parameter estimates at each iteration, as returned by
#' optim_procedure
#' @param vec_X Vectorised spot-level gene expression matrix
#' @param svd_obj Output of performing SVD on the spatial weight matrix. Also
#' includes slots for the inverses of U and V, as well as t(U) * V.
#' @param D_list List of expanded proportions of length K. Each is a diagonal
#' \eqn{np \times np} matrix whose values are the k'th cell type proportions,
#' repeated p times.
#' @param W Spatial weight matrix
#'
#' @return Vector of standard errors for each parameter
#'
#' @keywords internal

se_from_estim <- function(estim,
                          vec_X,
                          svd_obj,
                          D_list,
                          W) {
  K <- length(D_list)
  p <- (-3 + sqrt(9 + 8*ncol(estim)/K))/2

  params_list <- estim_to_cov(estim[nrow(estim),], K = K, p = p)
  P_list <- params_list[c(1:K)]
  Sigma_list <- params_list[c((K+1):(2*K))]

  Ainv_list <- get_ainv(P_list,
                        d = svd_obj$d,
                        U_inv = svd_obj$u_inv,
                        V_inv = svd_obj$v_inv,
                        tUxV = svd_obj$tUxV)
  V <- calcul_V_P(Ainv_list = Ainv_list,
                  Sigma_list = Sigma_list,
                  D_list = D_list)
  dl_dVinv <- calcul_dl_dVinv(vec_X = vec_X, V = V)
  dl_dVinv_Vinv <- t(solve(t(as.matrix(V)), t(as.matrix(dl_dVinv))))

  second_deriv_mat <- -calcul_d2l_dVinvdX_matrix(Ainv_list = Ainv_list,
                                                 D_list = D_list,
                                                 W = W,
                                                 Sigma_list = Sigma_list,
                                                 V = V,
                                                 dl_dVinv_Vinv = dl_dVinv_Vinv)
  hess_inv <- solve(second_deriv_mat)
  return(sqrt(diag(hess_inv)))
}
