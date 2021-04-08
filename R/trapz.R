#' Trapezoidal integration
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
trapz = function (x, y) {
  if(length(x) == 0)
    return(0)
  idx = 2:length(x)
  return(
    as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2
  )
}
