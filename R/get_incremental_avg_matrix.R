#' Calculate incremental average matrix
#'
#' @return incremental average matrix from development triangle
#' @param dTriangle development triangle
#' @export
get_incremental_avg_matrix <- function(dTriangle) {
  size <- nrow(dTriangle)
  cbind(dTriangle[, 1], (dTriangle[, (2:size)] + 0 *
                           dTriangle[, (1:(size - 1))]) -
          (dTriangle[, (1:(size - 1))] + 0 * dTriangle[, (2:size)]))
}
