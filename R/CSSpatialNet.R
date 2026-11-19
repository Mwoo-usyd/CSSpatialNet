#' CSSpatialNet
#'
#' @description
#' Description of function
#'
#' @param X Spot-level gene expression matrix, dimension \eqn{n \times p}, where
#' n is the number of spots and p is the number of genes
#' @param coords Spot-level spatial coordinates matrix, dimension
#' \eqn{n \times 2}
#' @param W Spatial weight matrix. If it is provided, it is used in the
#' optimisation procedure. If not, it is generated with the inputs coords and
#' lattice_shape
#' @param proportions Matrix of cell type proportions within each spot. Each row
#' represents one spot, and contains the proportions of that spot that is each
#' cell type. Dimension \eqn{n \times K}, where K is the number of cell types
#' @param tol Tolerance level for optimisation procedure. Optimisation steps
#' when absolute relative change in all parameters and log-likelihood from one
#' iteration to the next is within tol
#' @param max_iter Maximum number of iterations in optimisation procedure. Once
#' this is reached, the last estimates are returned, even if the absolute
#' relative changes of the parameters are not within tol
#' @param print_every Frequency of printing estimation procedure progress
#' @param alpha Alpha value in ADAM optimisation method
#' @param beta1 Beta1 value in ADAM optimisation method
#' @param beta2 Beta2 value in ADAM optimisation method
#' @param epsilon Epsilon value in ADAM optimisation method
#' @param lattice_shape String, either "hexagon" or "square", specifying if the
#' spots are on a hexagonal grid or a square grid.
#'
#' @return List of two lists, where the first list is titled estims and the
#' second list is titled SEs. The estims entry is a list of matrices of length
#' 2*K, where the first K matrices are the estimated cell-type-specific diagonal
#' matrices of spatial autocorrelation, and the last K matrices are the
#' estimated cell-type-specific gene-gene covariance matrices. The SEs slot
#' contains a list of 2*K matrices as well, where each entry represents the
#' standard error of that parameter. The structure is the same as the matrices
#' in the estims slot.
#'
#' @examples
#' data(SimulatedData)
#' CSSpatialNet_output <- CSSpatialNet(X = X,
#'   W = weights,
#'   proportions = props_mat)
#'
#' @export
#'



CSSpatialNet <- function(X,
                         coords,
                         W = NULL,
                         proportions,
                         tol = 0.01,
                         max_iter = 1000,
                         print_every = 200,
                         alpha = 0.01,
                         beta1 = 0.9,
                         beta2 = 0.999,
                         epsilon = 1e-8,
                         lattice_shape = "hexagon") {
  # create weights
  if (is.null(W)) {
    W <- create_weights(coords = coords, lattice_shape = lattice_shape)
  }
  svd_obj <- svd(W)
  svd_obj$v_inv <- solve(svd_obj$v)
  svd_obj$u_inv <- solve(svd_obj$u)
  svd_obj$tUxV <- t(svd_obj$u) %*% svd_obj$v

  prop_list <- asplit(proportions, MARGIN = 2)
  p <- ncol(X)
  D_list <- lapply(prop_list, function(x){rep(x,p)})
  D_list <- lapply(D_list, diag)

  # perform optimisation
  optim_out <- optim_procedure(W = W,
                               svd_obj = svd_obj,
                               vec_X = c(X),
                               tol = tol,
                               max_iter = max_iter,
                               print_every = print_every,
                               D_list = D_list,
                               alpha = alpha,
                               beta1 = beta1,
                               beta2 = beta2,
                               epsilon = epsilon)
  final_estims <- estim_to_cov(long_estim = optim_out[nrow(optim_out),],
                               K = ncol(proportions),
                               p = ncol(X))

  # calculate SEs
  SEs <- se_from_estim(estim = optim_out,
                       vec_X = c(X),
                       svd_obj = svd_obj,
                       D_list = D_list,
                       W = W)
  final_ses <- estim_to_cov(SEs,
                            K = ncol(proportions),
                            p = ncol(X))

  # return estimates and SEs
  return(list(estims = final_estims,
              SEs = final_ses))
}
