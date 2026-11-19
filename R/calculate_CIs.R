#' calculate_CIs
#'
#' @description
#' Function to calculate confidence intervals for each parameter in the
#' CSSpatialNet output, of a given confidence level.
#'
#' @param X Spot-level gene expression matrix, dimension \eqn{n \times p}, where
#' n is the number of spots and p is the number of genes
#' @param proportions Matrix of cell type proportions within each spot. Each row
#' represents one spot, and contains the proportions of that spot that is each
#' cell type. Dimension \eqn{n \times K}, where K is the number of cell types
#' @param CSSpatialNet_output Output from CSSpatialNet function performed on the
#' given X and proportions
#' @param conf_level Confidence level to be used when constructing the
#' confidence intervals. Default is 0.95
#'
#' @return Matrix of confidence intervals for each parameter. Each row
#' represents one parameter, and there are 7 columns, each as follows.
#' 1) Name of the cell type corresponding to that parameter
#' 2) Name of the gene (or pair of genes) corresponding to that parameter
#' 3) Parameter estimate from CSSpatialNet
#' 4) Standard error of parameter estimate from CSSpatialNet
#' 5) Parameter type, either spatial autocorrelation, variance, or covariance
#' 6) Lower bound of confidence interval
#' 7) Upper bound of confidence interval
#'
#' @importFrom vctrs vec_rep_each
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom stats qnorm
#' @importFrom rlang .data
#'
#' @examples
#' data(SimulatedData)
#' CSSpatialNet_output <- CSSpatialNet(X = X,
#'   W = weights,
#'   proportions = props_mat)
#' CIs_dataframe <- calculate_CIs(X = X,
#'   proportions = props_mat,
#'   CSSpatialNet_output = CSSpatialNet_output)
#'
#' @export

calculate_CIs <- function(X,
                          proportions,
                          CSSpatialNet_output,
                          conf_level = 0.95) {
  p <- ncol(X)
  if (is.null(colnames(X))) {
    gene_names <- paste0("Gene ", c(1:p))
  } else {
    gene_names <- colnames(X)
  }

  K <- ncol(proportions)
  if (is.null(colnames(proportions))) {
    cell_type_names <- paste0("Cell type ", c(1:K))
  } else {
    cell_type_names <- colnames(proportions)
  }

  long_cell_type_names <- c(vec_rep_each(cell_type_names, p),
                            vec_rep_each(cell_type_names, p*(p+1)/2))
  gene_pairs <- outer(gene_names, gene_names, function(a,b){paste0(a, ", ",b)})
  diag(gene_pairs) <- gene_names
  gene_pairs_ordered <- gene_pairs[lower.tri(gene_pairs, diag = TRUE)]
  long_gene_names <- c(rep(gene_names, K),
                       rep(gene_pairs_ordered, K))

  param_name_matrix <- matrix("Covariance", nrow = p, ncol = p)
  diag(param_name_matrix) <- rep("Variance", p)
  param_types_long <- c(rep("Spatial autocorrelation", p*K),
    rep(param_name_matrix[lower.tri(param_name_matrix, diag = TRUE)], K))

  long_estims <- matrix_to_estim(CSSpatialNet_output$estims)
  long_SEs <- matrix_to_estim(CSSpatialNet_output$SEs)

  CIs_df <- data.frame(cell_type = long_cell_type_names,
                       gene = long_gene_names,
                       estim = long_estims,
                       SE = long_SEs,
                       param_type = param_types_long)
  CIs_df <- CIs_df %>% mutate(lower = .data$estim +
                                qnorm((1-conf_level)/2)*.data$SE,
                              upper = .data$estim +
                                qnorm((1+conf_level)/2)*.data$SE)
  return(CIs_df)
}

