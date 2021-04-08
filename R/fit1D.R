#' 1D peak fit along CV axis
#'
#' @param mz0
#' @param CV0
#' @param dmz
#' @param dCV
#' @param mz
#' @param CV
#' @param MS
#' @param fwhm_cv_nom
#' @param weighted
#' @param refine_CV0
#' @param correct_overlap
#'
#' @return
#' @export
#'
#' @examples
fit1D <- function(
  mz0, CV0,
  dmz, dCV,
  mz, CV, MS,
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
  MSloc  = MS[selCV, selMz]
  mzloc  = mz[selMz]
  mMS    = apply(MSloc, 1, function(x) trapz(mzloc, x))

  mz_max = which.max(MSloc[which.max(mMS),])
  mz0 = mz[selMz][mz_max]
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz
  selMz  = mz >= mz1 & mz <= mz2 # Select mz area

  # Normal fit
  mzloc = mz[selMz]
  mMS = apply(MS[selCV, selMz], 1, function(x) trapz(mzloc,x))

  # Weighting scheme
  weights = rep(1,length(mMS))
  if(weighted){ # Poisson
    weights = 1 / mMS
    vm = min(mMS[mMS>0])
    weights[!is.finite(weights)] = 1 / vm
  }

  # First pass
  sigma0 = fwhm_cv_nom/2.355
  A0 = sqrt(2*pi) * sigma0 * max(mMS)
  lower = c(
    mu    = CV0 - dCV/10,
    sigma = 0.8 * sigma0,
    A     = 0.5 * A0
  )
  start = c(
    mu    = CV0,
    sigma = sigma0,
    A     = A0
  )
  upper = c(
    mu    = CV0 + dCV/10,
    sigma = 1.2 * sigma0,
    A     = 1.5 * A0
  )

  res = try(
    nls(
      mMS ~ A/(sqrt(2*pi)*sigma)*exp(-1/2*(CVf-mu)^2/sigma^2),
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
          mMS ~ A/(sqrt(2*pi)*sigma)*exp(-1/2*(CVf-mu)^2/sigma^2),
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
    # if(class(res) != 'try-error')
    #   print(cbind(lower,coef(summary(res))[,"Estimate"],upper))
  }

  return(
    list(
      mz0 = mz0,
      mz1 = mz1,
      mz2 = mz2,
      res = res,
      mMS = mMS,
      mMStot = mMStot,
      CVf = CVf
    )
  )
}
