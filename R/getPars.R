#' Get peak fit parameters
#'
#' @param res
#' @param dimfit
#'
#' @return
#' @export
#'
#' @examples
getPars = function(res, dimfit){
  if (dimfit == 2)
    getPars2D(res)
  else
    getPars1D(res, dimfit)
}
