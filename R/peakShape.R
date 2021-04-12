#' Calculate 1D Gaussian peak shape
#'
#' @param x (vector) a vector of coordinales
#' @param p (named vector) a named vector of parameters
#'
#' @return A vector of peak intensity values
#'
#' @export
#'
#' @examples
peakShape = function(
  x,
  p
) {
  if(length(p)==3)
    p['A'] / (sqrt(2*pi) * p['sigma'] ) *
    exp(-1/2*(x-p['mu'])^2/p['sigma']^2)
  else
    p['A'] / (sqrt(2*pi) * p['sy'] ) * # Marginalized on x
    exp(-1/2*(x-p['my'])^2/p['sy']^2)
}
