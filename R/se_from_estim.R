se_from_estim <- function(estim, z, svd_obj, D_list, W, new = TRUE) {
  K <- length(D_list)
  p <- (-3 + sqrt(9 + 8*ncol(estim)/K))/2 # rearrange formula number_of_params_to_estimate = K*p*(p+3)/2

  params_list <- estim_to_cov(estim[nrow(estim),], K = K, p = p)
  P_list <- params_list[c(1:K)]
  Sigma_list <- params_list[c((K+1):(2*K))]

  Ainv_list <- get_ainv(P_list,
                        d = svd_obj$d,
                        U_inv = svd_obj$u_inv,
                        V_inv = svd_obj$v_inv,
                        tUxV = svd_obj$tUxV)
  current_V <- calcul_V_P(current_Ainv_list = Ainv_list,
                          current_Sigma_list = Sigma_list,
                          D_list = D_list)
  dl_dVinv <- calcul_dl_dVinv(z = z, current_V = current_V)
  dl_dVinv_Vinv <- t(solve(t(current_V), t(dl_dVinv)))

  second_deriv_mat <- -calcul_d2l_dVinvdX_matrix(Ainv_list = Ainv_list,
                                                 D_list = D_list,
                                                 W = W,
                                                 Sigma_list = Sigma_list,
                                                 V = current_V,
                                                 dl_dVinv_Vinv = dl_dVinv_Vinv) # wrt neg log-likelihood
  hess_inv <- solve(second_deriv_mat)
  return(sqrt(diag(hess_inv)))
}
