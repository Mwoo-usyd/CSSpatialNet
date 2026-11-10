check_same_cell_type <- function(a,b,p,K) {
  # for entry a,b in matrix of second derivatives, does it correspond to two parameters of the same cell type (whether a Pk or Sigmak)?
  Pk_cell_types <- vec_rep_each(c(1:K), p)
  Sigmak_cell_types <- vec_rep_each(c(1:K), p*(p+1)/2)
  all_cell_types <- c(Pk_cell_types, Sigmak_cell_types)
  return(all_cell_types[a] == all_cell_types[b])
}
