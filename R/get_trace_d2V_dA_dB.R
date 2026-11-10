#' get_trace_d2V_dA_dB
#'
#' @description
#' Function to generate the component of the second derivative of the
#' log-likelihood that is only non-zero if the two parameters are of the same
#' cell type.
#'
#' @param a Index of first parameter
#' @param b Index of second parameter
#' @param p number of genes
#' @param K number of cell types
#' @param d2V_dPk_dPk_list List of second derivatives of V (covariance matrix
#' of vec(X)) with respect to two spatial autocorrelation parameters
#' @param d2V_dSigmak_dPk_list List of second derivatives of V (covariance
#' matrix of vec(X)) with respect to one spatial autocorrelation parameter and
#' one gene variance or covariance parameter
#' @param Vinv_dl_dVinv_Vinv V^{-1} times derivative of log-likelihood with
#' respect to V^{-1} times V^{-1}
#' @param all_cell_types vector of cell types (labelled 1 to K) for
#' each parameter
#'
#' @return Matrix of parts of second derivative of V that are non-zero only if
#' the two parameters are of the same cell type
#'
#' @keywords internal
#'

get_trace_d2V_dA_dB <- function(a,
                                b,
                                p,
                                K,
                                d2V_dPk_dPk_list,
                                d2V_dSigmak_dPk_list,
                                Vinv_dl_dVinv_Vinv,
                                all_cell_types) {
  if (check_same_cell_type(a = a, b = b, p = p, K = K)) {
    current_cell_type <- all_cell_types[a]
    param1 <- ifelse(a <= K*p, "Pk", "Sigmak")
    param2 <- ifelse(b <= K*p, "Pk", "Sigmak")
    a_adj <- ifelse(param1 == "Pk",
                    a - (current_cell_type-1)*p,
                    a - K*p - (current_cell_type-1)*p*(p+1)/2)
    b_adj <- ifelse(param2 == "Pk",
                    b - (current_cell_type-1)*p,
                    b - K*p - (current_cell_type-1)*p*(p+1)/2)
    #test_2nd_deriv_mat[a,b] <- paste0(a_adj, ", ", b_adj)

    if (param1 == "Pk" & param2 == "Pk") {
      lower <- min(a_adj,b_adj)
      upper <- max(a_adj,b_adj)
      list_position_d2V_dPk_dPk <- upper*(upper-1)/2 + lower + (current_cell_type-1) * p * (p+1)/2
      current_d2V_dPk_dPk <- d2V_dPk_dPk_list[[list_position_d2V_dPk_dPk]]
      trace_result <- calcul_trace_2input(a = Vinv_dl_dVinv_Vinv, b = current_d2V_dPk_dPk)
    } else if ((param1 == "Sigmak" & param2 == "Pk") | (param1 == "Pk" & param2 == "Sigmak")) {
      Pk_pos <- min(a,b)
      Sigmak_pos <- max(a,b)
      Pk_pos_adj <- Pk_pos - (current_cell_type - 1) * p
      Sigmak_pos_adj <- Sigmak_pos - K*p - (current_cell_type - 1) * p * (p+1)/2
      prev_cell_types_skip <- (current_cell_type-1) * p * p * (p+1)/2
      list_position_d2V_dSigmak_dPk <- prev_cell_types_skip + (Pk_pos_adj - 1) * p * (p+1)/2 + Sigmak_pos_adj
      current_d2V_dSigmak_dPk <- d2V_dSigmak_dPk_list[[list_position_d2V_dSigmak_dPk]]
      trace_result <- calcul_trace_2input(a = Vinv_dl_dVinv_Vinv, b = current_d2V_dSigmak_dPk)
    } else {trace_result <- 0}
  } else {trace_result <- 0}
  return(-trace_result)
}
