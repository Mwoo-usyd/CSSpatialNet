CSSpatialNet <- function(X,
                         coords,
                         proportions,
                         tol = 0.01,
                         max_iter = 1000,
                         print_every = 200,
                         alpha = 0.01,
                         beta1 = 0.9,
                         beta2 = 0.999,
                         epsilon = 1e-8) {
  # create weights
  weight_mat <- create_sparse_weights_grid(coords = coords, max_eigen_1 = TRUE, neighbours = "immediate")
  svd_obj <- svd(weight_mat)
  svd_obj$v_inv <- solve(svd_obj$v)
  svd_obj$u_inv <- solve(svd_obj$u)
  svd_obj$tUxV <- t(svd_obj$u) %*% svd_obj$v

  prop_list <- asplit(proportions, MARGIN = 2)
  D_list <- lapply(prop_list, function(x){rep(x,p)})
  D_list <- lapply(D_list, diag)

  # perform optimisation
  optim_out <- optim_procedure(W = weight_mat,
                               svd_obj = svd_obj,
                               x = c(X),
                               tol = tol,
                               max_iter = max_iter,
                               print_every - print_every,
                               D_list = D_list,
                               alpha = alpha,
                               beta1 = beta1,
                               beta2 = beta2,
                               epsilon = epsilon)

  # calculate SEs
  SEs <- se_from_estim(estim = optim_out, vec_X = c(X), svd_obj = svd_obj, D_list = D_list, W = weight_mat, new = TRUE) # need se_from_estim

  # return estimates and SEs
}
