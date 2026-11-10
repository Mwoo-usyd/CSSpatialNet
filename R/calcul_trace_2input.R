calcul_trace_2input <- function(a,b,dl_dVinv_Vinv = NULL) {
  if (is.null(dl_dVinv_Vinv)) {
    result <- a %*% b
  }
  else {
    result <- dl_dVinv_Vinv %*% a %*% b
  }
  return(sum(diag(result)))
}
