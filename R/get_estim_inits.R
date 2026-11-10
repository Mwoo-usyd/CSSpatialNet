#' get_estim_inits
#'
#' @description
#' Function to generate random initial parameters for CSSpatialNet optimisation
#' procedure
#'
#' @param p number of genes
#' @param K number of cell types
#'
#' @return vector of initial parameter values
#'
#' @keywords internal
#'

get_estim_inits <- function(p,K) {
  triangle_num <- p*(p-1)/2
  P_inits <- runif(K*p)
  Sigma_diag_inits <- runif(p, min = 1, max = 2)
  Sigma_offdiag_inits <- runif(triangle_num)
  Sigma_init_mat <- matrix(NA, nrow = p, ncol = p)
  Sigma_init_mat[lower.tri(Sigma_init_mat)] <- Sigma_offdiag_inits
  diag(Sigma_init_mat) <- Sigma_diag_inits
  estim_vals <- matrix(c(P_inits,
                         rep(Sigma_init_mat[lower.tri(Sigma_init_mat,
                                                      diag = TRUE)], K)),
                       nrow = 1, ncol = p*K*2+K*triangle_num)
  return(estim_vals)
}
