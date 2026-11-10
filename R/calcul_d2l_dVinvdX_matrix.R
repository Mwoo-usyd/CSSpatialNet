#' calcul_d2l_dVinvdX_matrix
#'
#' @description
#' This functions calculates the second derivative of the log-likelihood with
#' respect to V^{-1} (inverse of the covariance matrix of vec(X))
#' and each parameter.
#'
#' @param Ainv_list List of matrices of the form
#' $(I_{n} - P_{k,i} \otimes W)^{-1}$ where $P_{k,i}$ is the spatial autocorrelation
#' of the i'th gene in the k'th cell type
#' @param D_list List of expanded proportions of length K. Each is a diagonal
#' np *np matrix whose values are the k'th cell type proportions,
#' repeated p times.
#' @param W Spatial weight matrix
#' @param Sigma_list List of cell-type-specific covariance matrices
#' @param V Covariance matrix of vec(X)
#' @param dl_dVinv_Vinv Derivative of log-likelihood with respect to $V^{-1}$,
#' multiplied by $V^{-1}$
#'
#' @return Matrix of partial second derivatives of log-likelihood with respect
#' to $V^{-1}$ and each parameter.
#'
#' @keywords internal
#'
calcul_d2l_dVinvdX_matrix <- function(Ainv_list,
                                      D_list,
                                      W,
                                      Sigma_list,
                                      V,
                                      dl_dVinv_Vinv) {
  n <- nrow(W)
  K <- length(Ainv_list)
  p <- nrow(V)/n
  hessian_dimension <- K*p*(p+3)/2 # Sigmak and Pk matrices
  D_Ainv_list <- lapply(c(1:K), function(k){
    return(Ainv_list[[k]] * diag(D_list[[k]]))
  })
  Ainv_D_list <- lapply(c(1:K), function(k){
    return(t(t(Ainv_list[[k]]) * diag(D_list[[k]])))
  })

  # Sigmak matrices
  dV_dSigmak_list <- calcul_dV_dSigmak(Ainv_list = Ainv_list,
                                       D_list = D_list,
                                       n = n)
  t_dV_dSigmak_list <- lapply(dV_dSigmak_list, t)
  dV_dSigmak_x_Vinv_list <- lapply(t_dV_dSigmak_list, function(x){
    result <- t(solve(V, x))
    return(result)
  })
  t_dV_dSigmak_x_Vinv_list <- lapply(dV_dSigmak_list, function(x){
    result <- t(solve(V, x))
    return(result)
  })

  # Pk matrices
  dV_dPk_list <- calcul_dV_dPk(Ainv_list = Ainv_list,
                               D_list = D_list,
                               W = W,
                               Sigma_list = Sigma_list)
  t_dV_dPk_list <- lapply(dV_dPk_list, t)
  dV_dPk_x_Vinv_list <- lapply(t_dV_dPk_list, function(x){
    result <- t(solve(V,x))
    return(result)
  })
  t_dV_dPk_x_Vinv_list <- lapply(dV_dPk_list, function(x){
    result <- t(solve(V,x))
    return(result)
  })

  dV_dX_Vinv_list <- append(dV_dPk_x_Vinv_list, dV_dSigmak_x_Vinv_list)
  t_dV_dX_Vinv_list <- append(t_dV_dPk_x_Vinv_list, t_dV_dSigmak_x_Vinv_list)

  # extract traces for first line
  grid_params <- expand.grid(i = seq_len(hessian_dimension), j = seq_len(hessian_dimension))
  first_line_flat_results <- mapply(function(i, j) calcul_trace_2input(a = t_dV_dX_Vinv_list[[i]],
                                                                       b = dV_dX_Vinv_list[[j]]),
                                    grid_params$i, grid_params$j, SIMPLIFY = FALSE)
  first_line_result_matrix <- matrix(as.numeric(first_line_flat_results), nrow = hessian_dimension, ncol = hessian_dimension) * (-0.5)

  # cross line 1
  cross1_flat_results <- mapply(function(i, j) calcul_trace_2input(a = dV_dX_Vinv_list[[i]],
                                                                   b = dV_dX_Vinv_list[[j]],
                                                                   dl_dVinv_Vinv = dl_dVinv_Vinv),
                                grid_params$i, grid_params$j, SIMPLIFY = FALSE)
  cross1_result_matrix <- matrix(as.numeric(cross1_flat_results), nrow = hessian_dimension, ncol = hessian_dimension)
  # cross line 2
  cross2_flat_results <- mapply(function(i, j) calcul_trace_2input(a = dV_dX_Vinv_list[[j]],
                                                                   b = dV_dX_Vinv_list[[i]],
                                                                   dl_dVinv_Vinv = dl_dVinv_Vinv),
                                grid_params$i, grid_params$j, SIMPLIFY = FALSE)
  cross2_result_matrix <- matrix(as.numeric(cross2_flat_results), nrow = hessian_dimension, ncol = hessian_dimension)

  # within-cell-type derivatives

  within_deriv_mat <- matrix(0, nrow = hessian_dimension, ncol = hessian_dimension)
  Pk_cell_types <- vec_rep_each(c(1:K), p)
  Sigmak_cell_types <- vec_rep_each(c(1:K), p*(p+1)/2)
  all_cell_types <- c(Pk_cell_types, Sigmak_cell_types)
  Vinv_dl_dVinv_Vinv <- solve(V, dl_dVinv_Vinv)

  # create d2V_dPk_dPk
  d2V_dPk_dPk_list <- calcul_d2V_dPk_dPk(Ainv_list = Ainv_list, D_Ainv_list = D_Ainv_list, Ainv_D_list = Ainv_D_list, W = W, Sigma_list = Sigma_list)
  # create d2v_dPk_dSigmak
  d2V_dSigmak_dPk_list <- calcul_d2V_dSigmak_dPk(Ainv_list = Ainv_list, D_list = D_list, W = W, D_Ainv_list = D_Ainv_list, Ainv_D_list = Ainv_D_list)

  # fill matrix of second derivatives specifically for the within-cell-type term
  flat_results_d2V_dA_dB <- mapply(function(i, j) get_trace_d2V_dA_dB(a = i,
                                                                      b = j,
                                                                      p = p,
                                                                      K = K,
                                                                      d2V_dPk_dPk_list = d2V_dPk_dPk_list,
                                                                      d2V_dSigmak_dPk_list = d2V_dSigmak_dPk_list,
                                                                      Vinv_dl_dVinv_Vinv = Vinv_dl_dVinv_Vinv,
                                                                      all_cell_types = all_cell_types),
                                   grid_params$i, grid_params$j, SIMPLIFY = FALSE)
  matrix_results_d2V_dA_dB <- matrix(as.numeric(flat_results_d2V_dA_dB), nrow = hessian_dimension, ncol = hessian_dimension)

  # final result
  return(first_line_result_matrix + cross1_result_matrix + cross2_result_matrix + matrix_results_d2V_dA_dB)
}
