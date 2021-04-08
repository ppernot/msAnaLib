#' 1D peak fit along m/z axis
#'
#' @param mz0
#' @param CV0
#' @param dmz
#' @param dCV
#' @param mz
#' @param CV
#' @param MS
#' @param fwhm_mz_nom
#' @param weighted
#'
#' @return
#' @export
#'
#' @examples
fit1D_MS <- function(
  mz0, CV0,
  dmz, dCV,
  mz, CV, MS,
  fwhm_mz_nom,
  weighted = FALSE
) {

  # Select CV0
  iCV = which(CV >= CV0)[1]

  # Select mz window
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz
  selMz  = mz >= mz1 & mz <= mz2 # Select mz area
  # mz profile
  MStot = MS[iCV,]
  mzl = as.vector(mz[selMz])
  MSl = MS[iCV, selMz]

  # Normal fit
  weights = rep(1,length(MSl))
  if(weighted){ # Poisson
    weights = 1 / MSl
    vm = min(MSl[MSl>0])
    weights[!is.finite(weights)] = 1 / vm
  }

  # First pass
  mu = mzl[which.max(MSl)]
  sigma0 = fwhm_mz_nom/2.355
  A0 = sqrt(2*pi) * sigma0 * max(MSl)
  lower = c(
    mu    = mu - dmz,
    sigma = 0.1 * sigma0,
    A     = 0.5 * A0
  )
  start = c(
    mu    = mu,
    sigma = sigma0,
    A     = A0
  )
  upper = c(
    mu    = mu + dmz,
    sigma = 2.0 * sigma0,
    A     = 1.5 * A0
  )

  res = try(
    nls(
      MSl ~ A/(sqrt(2*pi)*sigma)*exp(-1/2*(mzl-mu)^2/sigma^2),
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
          MSl ~ A/(sqrt(2*pi)*sigma)*exp(-1/2*(mzl-mu)^2/sigma^2),
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
      mzl = mzl,
      MSl = MSl,
      res = res,
      mz0 = mz0,
      mz1 = mz1,
      mz2 = mz2,
      CVf = CV,
      mMStot = MStot
    )
  )
}
