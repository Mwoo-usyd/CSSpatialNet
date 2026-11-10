get_estim_inits <- function(p,K) {
  triangle_num <- p*(p-1)/2
  P_inits <- runif(K*p)
  Sigma_diag_inits <- runif(p, min = 1, max = 2)
  Sigma_offdiag_inits <- runif(triangle_num)
  Sigma_init_mat <- matrix(NA, nrow = p, ncol = p)
  Sigma_init_mat[lower.tri(Sigma_init_mat)] <- Sigma_offdiag_inits
  # Sigma_init_mat[upper.tri(Sigma_init_mat)] <- Sigma_offdiag_inits # don't need if just taking lower.tri(diag = TRUE) anyway?
  diag(Sigma_init_mat) <- Sigma_diag_inits
  estim_vals <- matrix(c(P_inits,
                         rep(Sigma_init_mat[lower.tri(Sigma_init_mat, diag = TRUE)], K)),
                       nrow = 1, ncol = p*K*2+K*triangle_num)
  return(estim_vals)
}
