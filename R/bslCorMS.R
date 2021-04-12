#' Correct baseline
#'
#' @param MS (matrix; mandatory) a 2D matrix
#' @param baseline_cor (string; mandatory) the correction method to apply
#'   (default: "median+global"). At the moment, the options are:
#'   \describe{
#'     \item{'median' or 'mean'}{The constant baseline correction is
#'       based on the median (default) or mean of the signal.}
#'     \item{'row' or 'column' or 'global'}{The correction is estimated
#'       from and applied  to each row or column of the matrix, or
#'       estimated from and applied to the whole matrix (default).}
#'     \item{NULL}{The baseline correction is skipped.}
#'   }
#'
#' @return A baseline-corrected matrix.
#' @export
#'
#' @examples
bslCorMS = function(MS, baseline_cor = 'median+global') {
  # Baseline correction in the hypothesis of a constant shift.
  # If a constant seems insufficient, consider using
  # one of the methods in "baseline" package...

  if(!is.null(baseline_cor)) {

    if(grepl('median',baseline_cor)) {
      cor_fun = median
    } else if(grepl('mean',baseline_cor)) {
      cor_fun = mean
    } else {
      message('>>> Error in bslCorMS:\n>>>>>> Unknown method: ',baseline_cor)
      stop(call.=FALSE)
    }

    if (grepl('row', baseline_cor)) {
      cat('>>> Row-wise baseline correction...\n')
      for (i in 1:nrow(MS))
        MS[i, ] = MS[i, ] - cor_fun(MS[i, ])
    } else if (grepl('col', baseline_cor)) {
      cat('>>> Column-wise baseline correction...\n')
      for (i in 1:ncol(MS))
        MS[, i] = MS[, i] - cor_fun(MS[, i])
    } else if (grepl('glob', baseline_cor)) {
      cat('>>> Global baseline correction...\n')
      MS = MS - cor_fun(MS)
    } else {
      message('>>> Error in bslCorMS:\n>>>>>> Unknown method: ',baseline_cor)
      stop(call.=FALSE)
    }

  }

  return(MS)
}
