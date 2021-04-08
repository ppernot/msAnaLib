#' Correct baseline
#'
#' @param MS
#' @param baseline_cor
#'
#' @return
#' @export
#'
#' @examples
bslCorMS = function(MS, baseline_cor = 'median') {
  # Baseline correction based on mode of the data distribution
  # in the hypothesis of a constant shift
  # FTICR has gamma noise distribution (doi:10.1016/j.aca.2009.10.043)
  # The median seems to be a good (constant) baseline estimator

  if(!is.null(baseline_cor)) {

    if(grepl('median',baseline_cor)) {
      cor_fun = median
    } else {
      cor_fun = mean
    }

    if(grepl('perCV',baseline_cor)) {
      cat('>>> Baseline correction of individual MSs...\n')
      # Operate correction on individual MSs
      for(i in 1:nrow(MS))
        MS[i,] = MS[i,]- cor_fun(MS[i,])
    } else {
      cat('>>> Global baseline correction...\n')
      # Global correction
      MS = MS - cor_fun(MS)
    }

  }

  return(MS)
}
