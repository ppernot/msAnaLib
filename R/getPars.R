#' Get peak fit parameters
#'
#' @param res (nls-object) a nls fit result
#' @param dimfit (integer) dimension of the fit
#'
#' @return A list of best-fit parameters:
#' \describe{
#'   \item{v}{vector of Gaussian peak best-fit parameters}
#'   \item{u_v}{uncertainty on v elements}
#'   \item{mzopt}{peak center along m/z}
#'   \item{u_mz}{uncertainty on mzopt}
#'   \item{cvopt}{peak center along CV}
#'   \item{u_cv}{uncertainty on cvopt}
#'   \item{fwhm_mz}{peak FWHM along m/z}
#'   \item{u_fwhm_mz}{uncertainty on fwhm_mz}
#'   \item{fwhm_cv}{peak FWHM along CV}
#'   \item{u_fwhm_cv}{uncertainty on fwhm_cv}
#'   \item{area}{peak area}
#'   \item{u_area}{uncertainty on area}
#' }
#' Depending on 'dimfit', some values might be NAs.
#'
#' @export
#'
#' @examples
getPars = function(
  res,
  dimfit = c(0,1,2)
){
  dimfit = match.arg(dimfit)
  if (dimfit == 2)
    getPars2D(res)
  else
    getPars1D(res, dimfit)
}
