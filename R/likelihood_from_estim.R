likelihood_from_estim <- function(estim, svd_obj, D_list, z, new = FALSE) {
  K <- length(D_list)
  p <- (-3 + sqrt(9 + 8*ncol(estim)/K))/2 # rearrange formula number_of_params_to_estimate = K*p*(p+3)/2
  n <- length(svd_obj$d)
  if (new) {
    params_list <- estim_to_cov_new(estim[nrow(estim),], K = K, p = p)
  } else {
    params_list <- estim_to_cov(estim[nrow(estim),], K = K, p = p)
  }

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
  likelihood_val <- dmvnorm(x = z, mean = rep(0, n*p), sigma = as.matrix(current_V), log = TRUE)
  return(likelihood_val)
}
