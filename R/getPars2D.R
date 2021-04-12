#' Extract parameters from 2D Gaussian peak fit
#'
#' @param res (nls-object) a nls fit result
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
#'
#' @export
#'
#' @examples
getPars2D <- function(res) {
  # Get best params and uncertainty
  v   = summary(res)$parameters[,"Estimate"]   # Best params
  u_v = summary(res)$parameters[,"Std. Error"] # Uncertainty

  # Transform params to quantities of interest
  mu1     = v['mx']
  u_mu1   = u_v['mx']
  mu2     = v['my']
  u_mu2   = u_v['my']
  fwhm1   = 2.355 * v['sx']
  u_fwhm1 = 2.355 * u_v['sx']
  fwhm2   = 2.355 * v['sy']
  u_fwhm2 = 2.355 * u_v['sy']
  area    = v['A']
  u_area  = u_v['A']

  return(
    list(
      v         = v,
      u_v       = u_v,
      mzopt     = mu1,
      u_mz      = u_mu1,
      cvopt     = mu2,
      u_cv      = u_mu2,
      fwhm_mz   = fwhm1,
      u_fwhm_mz = u_fwhm1,
      fwhm_cv   = fwhm2,
      u_fwhm_cv = u_fwhm2,
      area      = area,
      u_area    = u_area
    )
  )

}
