#' Get apparatus-based specifications for peaks
#'
#' @param ms_type
#'
#' @return
#' @export
#'
#' @examples
getPeakSpecs = function(ms_type) {

  if(! ms_type %in% c('esquire','fticr'))
    stop('>>> Bad ms_type !')

  if(ms_type == 'fticr') {

    # Acceptable range for MS peak width
    fwhm_mz_min = 0.001
    fwhm_mz_max = 0.02

    # Nominal value for peak fit
    fwhm_mz_nom = 0.01

    # Width of mz window around exact mz for signal averaging
    dmz = 5 * fwhm_mz_nom # ~ 10-sigma interval

    # Constant baseline over all matrix; median estimator
    baseline_cor = 'global+median'

  } else { # Default: Esquire

    # Acceptable range for MS peak width
    fwhm_mz_min = 0.1
    fwhm_mz_max = 0.5

    # Nominal value for peak fit
    fwhm_mz_nom = 0.2

    # Width of mz window around exact mz for signal averaging
    dmz = 5 * fwhm_mz_nom

    # Constant baseline over all matrix; median estimator
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
