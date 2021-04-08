#' Weighted mean and uncertainty
#'
#' @param x
#' @param ux
#'
#' @return
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
