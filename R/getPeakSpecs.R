#' Get apparatus-based specifications for peaks
#'
#' @param ms_type (string; optional) Type of mass spectrum,
#'   one of `esquire` (default) or `fticr`.
#'
#' @return A list with the following arguments:
#' \describe{
#'   \item{fwhm_mz_min, fwhm_mz_max}{Acceptable range for m/z peak width.}
#'   \item{fwhm_mz_nom}{Nominal value of m/z peak width.}
#'   \item{dmz}{Width of mz window around exact mz for signal averaging.}
#'   \item{baseline_cor}{Baseline correction method (see\link{bslCorMS})}
#'   \item{fwhm_cv_min, fwhm_cv_max}{Acceptable range for CV peak width.}
#'   \item{fwhm_cv_nom}{Nominal value of CV peak width.}
#'   \item{dCV}{Width of CV window around CV_ref for peak fit.}
#'   \item{area_min}{Minimal peak area for filtering.}
#' }
#'
#' @export
#'
#' @examples
getPeakSpecs = function(
  ms_type = c('esquire','fticr')
) {

  ms_type = match.arg(ms_type)

  if(ms_type == 'fticr') {

    # Acceptable range for MS peak width
    fwhm_mz_min = 0.001
    fwhm_mz_max = 0.02

    # Nominal value for peak fit
    fwhm_mz_nom = 0.01

    # Width of mz window around exact mz for signal averaging
    dmz = 5 * fwhm_mz_nom # ~ 10-sigma interval

    # Constant baseline over all matrix; median estimator seems
    # OK for the TI-CR spectra tested so far...
    baseline_cor = 'median+global'

  } else { # Default: Esquire

    # Acceptable range for MS peak width
    fwhm_mz_min = 0.1
    fwhm_mz_max = 0.5

    # Nominal value for peak fit
    fwhm_mz_nom = 0.2

    # Width of mz window around exact mz for signal averaging
    dmz = 5 * fwhm_mz_nom

    # Esquire MSs require no baseline correction
    baseline_cor = NULL

  }

  # Acceptable range for CV peak width
  fwhm_cv_min = 0.5
  fwhm_cv_max = 1.5

  # Nominal value for CV peak's FWHM constraint
  fwhm_cv_nom = 1.0

  # Width of CV window around reference CV for peak fit
  dCV = 2.5

  # Min area for validation of  peak
  area_min = 10

  return(
    list(
      fwhm_mz_min  = fwhm_mz_min,
      fwhm_mz_max  = fwhm_mz_max,
      fwhm_mz_nom  = fwhm_mz_nom,
      dmz          = dmz,
      baseline_cor = baseline_cor,
      fwhm_cv_min  = fwhm_cv_min,
      fwhm_cv_max  = fwhm_cv_max,
      fwhm_cv_nom  = fwhm_cv_nom,
      dCV          = dCV,
      area_min     = area_min
    )
  )

}
