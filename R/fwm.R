#' Weighted mean and uncertainty by the Dersimonian-Laird method.
#'
#' A dark uncertainty component is added for the statistical
#' compatibility of the averaged values.
#'
#' @param x (vector) set of values to be averages.
#' @param ux (vector) uncertainties of values to be averaged.
#'
#' @return A list with elements:
#' \describe{
#'   \item{wm}{Weighted mean.}
#'   \item{uwm}{Uncertainty of the weighted mean.}
#' }
#'
#' @export
#'
#' @examples
fwm = function(x, ux) {
  # Weigthed mean and its uncertainty using Cochran's ANOVA for weights
  # Ref: C. Rivier et al. (2014) Accredit. Qual. Assur. 19, 269–274
  # (doi:https://doi.org/10.1007/s00769-014-1066-3)

  # u2 is the part of the variance not explained by ux (parametric unc.)
  u2 = max(0, var(x) - mean(ux ^ 2)) # Cochran's ANOVA
  ut2 = u2 + ux^2                    # Combined variance for each point

  # Variance weights
  w = (1 / ut2) / sum(1 / ut2)

  # Standard weighted mean and uncertainty
  wm  = weighted.mean(x, w)
  uwm = sqrt(1 / sum(1 / ut2))

  return(list(wm = wm, uwm = uwm))
}
