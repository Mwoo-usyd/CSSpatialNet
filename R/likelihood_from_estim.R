#' likelihood_from_estim
#'
#' @description
#' Function to calculate the log-likelihood of the observed data under the
#' covariance matrix created with current parameter values
#'
#' @param estim Matrix of parameter estimates over all iterations as outputted
#' by optim_procedure
#' @param svd_obj Output of performing SVD on the spatial weight matrix. Also
#' includes slots for the inverses of U and V, as well as t(U) * V.
#'

likelihood_from_estim <- function(estim,
                                  svd_obj,
                                  D_list,
                                  vec_X) {
  K <- length(D_list)
  p <- (-3 + sqrt(9 + 8*ncol(estim)/K))/2
  n <- length(svd_obj$d)
  params_list <- estim_to_cov(estim[nrow(estim),], K = K, p = p)

  P_list <- params_list[c(1:K)]
  Sigma_list <- params_list[c((K+1):(2*K))]

  Ainv_list <- get_ainv(P_list,
                        d = svd_obj$d,
                        U_inv = svd_obj$u_inv,
                        V_inv = svd_obj$v_inv,
                        tUxV = svd_obj$tUxV)
  V <- calcul_V_P(current_Ainv_list = Ainv_list,
                  current_Sigma_list = Sigma_list,
                  D_list = D_list)
  likelihood_val <- dmvnorm(x = vec_X, mean = rep(0, n*p), sigma = as.matrix(V), log = TRUE)
  return(likelihood_val)
}
