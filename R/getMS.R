#' Get MS from file
#'
#' @param file (string; mandatory) Path to DMS file
#' @param ms_type (string; optional) Type of mass spectrum,
#'   one of `esquire` (default) or `fticr`.
#'
#' @return List of three elements: `mz` = vector of m/z coordinates,
#'   `time` = vector of time coordinates, and `MS` = DMS matrix.
#'
#' @export
#'
#' @examples
getMS = function(
  file,
  ms_type = c('esquire','fticr')
) {

  ms_type = match.arg(ms_type)

  cat('\n>>> Reading',ms_type,'MS in file:',file,'\n')
  MS0 = as.data.frame(
    data.table::fread( # Much faster than read.table
      file = file ,
      header = FALSE ,
      sep = ',' ,
      stringsAsFactors = FALSE
    )
  )
  time = MS0[, 1]

  if(ms_type == 'esquire') {

    range_mz = as.numeric(unlist(strsplit(MS0[1, 7], split = '-')))
    n_del_mz = as.numeric(unlist(strsplit(MS0[1, 8], split = '/')))
    nchan    = n_del_mz[1]
    del_mz   = n_del_mz[2]
    mz       = range_mz[1] + (0:(nchan - 1)) * del_mz
    MS       = as.matrix(MS0[, -(1:8)],
                         ncol = length(mz),
                         byrow = FALSE)

  } else {

    mySplit = function(x) as.numeric(unlist(strsplit(x, split = ' ')))
    mz = c()
    msLen = ncol(MS0)-8
    MS = matrix(NA, nrow = nrow(MS0), ncol = msLen)
    for (j in 1:nrow(MS0)) {
      dbl = vapply(
        MS0[j, 9:(msLen + 8)],
        FUN = mySplit,
        FUN.VALUE = numeric(2)
      )
      mz      = dbl[1, ]
      MS[j, ] = dbl[2, ]
    }

  }

  return(
    list(
      mz = mz,
      time = time,
      MS = MS
    )
  )
}
