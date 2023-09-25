#' 2D peak fit
#'
#' @param mz0 (numeric) Center of the m/z fit window.
#' @param CV0 (numeric) Center of the CV fit window
#' @param dmz (numeric) Width of the m/z fit window.
#' @param dCV (numeric) Width of the CV fit window.
#' @param mz (vector) m/z vector
#' @param CV (vector) CV vector
#' @param MS (matrix) DMS matrix
#' @param fwhm_mz_nom (numeric) Initial value for the m/z peak FWHM.
#' @param fwhm_cv_nom (numeric) Initial value for the CV peak FWHM.
#' @param weighted (logical; optional) Poisson weighting of data
#'   (default: FALSE).
#' @param refine_CV0 (logical; optional) Center the CV fit window on the
#'   max of the signal (default: FALSE).
#' @param correct_overlap (logical; optional) Use a correction
#'   for the peaks overlap (logical: FALSE). Experimental, to be
#'   implemented.
#'
#'
#' @return A list with elements:
#' \describe{
#'  \item{res}{a \link{nls} object containing the fit results}
#'  \item{mz0}{the center of the m/z fit window.}
#'  \item{mz1}{the lower bound of the m/z fit window.}
#'  \item{mz2}{the upper bound of the m/z fit window.}
#'  \item{CVf}{CV vector in fit window.}
#'  \item{mMS}{Averaged spectrum along the CV axis.}
#'  \item{mMStot}{Same as mMS ???}
#'  \item{mzloc}{m/z vector in the fit window.}
#'  \item{MSloc}{DMS signal in the fit window.}
#' }
#' @export
#'
#' @examples
fit2D <- function(
  mz0, CV0,
  dmz, dCV,
  mz, CV, MS,
  fwhm_mz_nom,
  fwhm_cv_nom,
  weighted = FALSE,
  refine_CV0 = FALSE,
  correct_overlap = FALSE
) {

  # Select mz window
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz
  selMz  = mz >= mz1 & mz <= mz2 # Select mz area
  # Integrated CV profile
  mzloc = mz[selMz]
  mMStot = apply(MS[, selMz], 1, function(x) trapz(mzloc,x))

  # Select CV window
  if(is.na(CV0))
    CV0 = CV[which.max(mMStot)]
  CV1 = CV0 - dCV/2
  CV2 = CV0 + dCV/2
  selCV = CV >= CV1 & CV <= CV2

  if(refine_CV0){
    im = which.max(mMStot[selCV])
    CV0 = CV[selCV][im]
    CV1 = CV0 - dCV/2
    CV2 = CV0 + dCV/2
    selCV = CV >= CV1 & CV <= CV2
  }
  CVf = CV[selCV]

  # Refine mz window
  MSloc = MS[selCV, selMz]
  mzloc = mz[selMz]
  mMS = apply(MSloc, 1, function(x) trapz(mzloc,x))
  mz_max = which.max(MSloc[which.max(mMS),])
  mz0 = mz[selMz][mz_max]
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz
  selMz  = mz >= mz1 & mz <= mz2 # Select mz area

  MSloc = MS[selCV, selMz]
  mzloc = mz[selMz]
  mMS = apply(MSloc, 1, function(x) trapz(mzloc,x))

  # Normal fit
  grid = expand.grid(x = mzloc, y = CVf)
  x = as.vector(grid$x)
  y = as.vector(grid$y)
  z = as.vector(t(MSloc))

  weights = rep(1,length(z))
  if(weighted){ # Poisson
    weights = 1 / z
    vm = min(z[z>0])
    weights[!is.finite(weights)] = 1 / vm
  }

  # First pass
  maxz = which.max(z)
  mx = x[maxz]
  sx0 = fwhm_mz_nom / 2.355
  sy0 = fwhm_cv_nom / 2.355
  A0  = 2 * pi * sx0 * sy0 * max(z)
  lower = c(
    mx  = mx - dmz / 2,
    sx  = sx0 / 4,
    my  = CV0 - dCV / 10,
    sy  = 0.8 * sy0,
    A   = 0.5 * A0
  )
  start = c(
    mx  = mx,
    sx  = sx0,
    my  = CV0,
    sy  = sy0,
    A   = A0
  )
  upper = c(
    mx  = mx + dmz / 2,
    sx  = 4 * sx0,
    my  = CV0 + dCV / 10,
    sy  = 1.5 * sy0,
    A   = 1.5 * A0
  )

  res = try(
    nls(
      z ~ A/(2*pi*sx*sy)*exp(
        -0.5 * ( (x-mx)^2/sx^2 + (y-my)^2/sy^2 ) ),
      start = start,
      lower = lower,
      upper = upper,
      algorithm = 'port',
      weights = weights,
      control = list(tol=1e-5,warnOnly=TRUE)
    ),
    silent = TRUE
  )

  if(class(res) != 'try-error') {
    if(res$convergence != 0) {
      # Attempt second pass from pass1 optimum
      start = as.list(coef(summary(res))[,"Estimate"])
      res = try(
        nls(
          z ~ A/(2*pi*sx*sy)*exp(
            -0.5 * ( (x-mx)^2/sx^2 + (y-my)^2/sy^2 ) ),
          start = start,
          lower = lower,
          upper = upper,
          algorithm = 'port',
          weights = weights,
          control = list(tol=1e-5)
        ),
        silent = TRUE
      )
    }
  }

  return(
    list(
      mz0 = mz0,
      mz1 = mz1,
      mz2 = mz2,
      res = res,
      mMS = mMS,
      mMStot = mMStot,
      CVf = CVf,
      mzloc = mzloc,
      MSloc = MSloc
    )
  )
}
