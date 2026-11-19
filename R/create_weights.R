#' create_weights
#'
#' @description
#' Function to create spatial weights from two-dimensional coordinates
#'
#' @param coords Matrix of spatial coordinates. Dimension \eqn{n \times 2}
#' @param lattice_shape String, either "hexagon" or "square", specifying if the
#' spots are on a hexagonal grid or a square grid.
#'
#' @return Spatial weight matrix of dimension n x n to be used for parameter
#' estimation.
#'
#' @importFrom stats dist
#'
#' @keywords internal

create_weights <- function(coords, lattice_shape = "hexagon") {
  n <- nrow(coords)
  distance_matrix <- as.matrix(dist(coords))
  smallest_distance <- min(distance_matrix[distance_matrix != 0])
  if (lattice_shape == "square") {
    radius <- smallest_distance * sqrt(2) + 10^-3
  } else if (lattice_shape == "hexagon") {
    radius <- smallest_distance + 10^-3
  } else {return("Invalid lattice shape")}
  weights <- matrix(0, nrow = n, ncol = n)
  weights[distance_matrix < radius] <- 1
  diag(weights) <- 0
  weights <- weights/rowSums(weights)
  weights <- 0.5 * (weights + t(weights))
  weights <- weights/max(abs(eigen(weights)$values))
  return(weights)
}
