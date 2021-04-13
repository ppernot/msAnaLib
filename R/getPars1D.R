#' Extract parameters from 1D peak fit
#'
#' @param res (nls-object) a nls fit result
#' @param fit_dim (integer) dimension of the fit (0 or 1)

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
getPars1D <- function(
  res,
  fit_dim = c(1,0)
) {

  if(!fit_dim %in% 0:1) {
    message('>>> Error: fit_dim should be 0 or 1 !')
    stop(call. = FALSE)
  }

  # Get best params and uncertainty
  v   = summary(res)$parameters[,"Estimate"]   # Best params
  u_v = summary(res)$parameters[,"Std. Error"] # Uncertainty

  # Transform params to quantities of interest
  mu     = v['mu']
  u_mu   = u_v['mu']
  fwhm   = 2.355 * v['sigma']
  u_fwhm = 2.355 * u_v['sigma']
  area   = v['A']
  u_area = u_v['A']

  if (fit_dim == 1)
    return(
      list(
        v         = v,
        u_v       = u_v,
        mzopt     = NA,
        u_mz      = NA,
        cvopt     = mu,
        u_cv      = u_mu,
        fwhm_mz   = NA,
        u_fwhm_mz = NA,
        fwhm_cv   = fwhm,
        u_fwhm_cv = u_fwhm,
        area      = area,
        u_area    = u_area
      )
    )
  else
    return(
      list(
        v         = v,
        u_v       = u_v,
        mzopt     = mu,
        u_mz      = u_mu,
        cvopt     = NA,
        u_cv      = NA,
        fwhm_mz   = fwhm,
        u_fwhm_mz = u_fwhm,
        fwhm_cv   = NA,
        u_fwhm_cv = NA,
        area      = area,
        u_area    = u_area
      )
    )

}
