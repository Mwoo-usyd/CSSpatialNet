estim_to_cov <- function(long_estim,
                         K,
                         p) {
  triangle_num <- p*(p-1)/2
  large_triangle <- (p+1)*p/2
  P_list <- list()
  for (current_k in 1:K) {
    vals_to_take <- c(((current_k-1)*p + 1):(current_k*p))
    P_list <- P_list %>% append(list(diag(long_estim[vals_to_take])))
  }
  # Extract Sigma_k matrices
  Sigma_list <- list()
  for (current_k in 1:K) {
    vals_to_take <- c(((current_k-1)*large_triangle + 1):(current_k*large_triangle)) + p*K
    current_Sigmak_entries <- long_estim[vals_to_take]
    current_Sigmak <- matrix(0, nrow = p, ncol = p)
    current_Sigmak[lower.tri(current_Sigmak, diag = TRUE)] <- current_Sigmak_entries # need to set upper triangle as well
    diag_entries <- diag(current_Sigmak)
    current_Sigmak <- current_Sigmak + t(current_Sigmak) - diag(diag_entries)
    Sigma_list <- Sigma_list %>% append(list(current_Sigmak))
  }
  return(c(P_list, Sigma_list))
}
