#' calcul_trace_2input
#'
#' @description
#' This function calculates the trace of the product of two matrices, since this
#' often needed in the second derivative calculations. There is the option to
#' include dl_dVinv_Vinv since the method often requires the trace of
#' dl_dVinv_Vinv multiplied with two other matrices.
#'
#' @param a The first matrix
#' @param b The second matrix
#' @param dl_dVinv_Vinv (Optional) the dl_dVinv_Vinv term in case you want the
#' trace of dl_dVinv_Vinv times a times b (matrix multiplication)
#'
#' @return Trace of a times b, or dl_dVinv_Vinv times a times b if dl_dVinv_Vinv
#' is specified
#'
#' @keywords internal

calcul_trace_2input <- function(a,
                                b,
                                dl_dVinv_Vinv = NULL) {
  if (is.null(dl_dVinv_Vinv)) {
    result <- a %*% b
  }
  else {
    result <- dl_dVinv_Vinv %*% a %*% b
  }
  return(sum(diag(result)))
}
