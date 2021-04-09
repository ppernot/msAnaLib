#' Down-sample a FT-ICR MS by choosing points closest to a preset
#' regular grid or around a set of targets.
#'
#' @param file_in (string; mandatory) Input DMS file to be compressed
#' @param file_out (string; mandatory) Output (compressed) file
#' @param ms_type (string; optional) Type of mass spectrum,
#'   one of `esquire` or `fticr` (default).
#'   Compression only for `fticr` MS.
#' @param compMode (string; optional) Compression mode:
#'   either `grid`, `targets` (default),
#'   or their combination `grid+targets`.
#' @param mzMin (numeric; optional) Min m/z value of compression grid
#' @param mzMax (numeric; optional) Max m/z value of compression grid
#' @param dmz (numeric; optional) m/z step or compression grid
#' @param targets (string; optional)  Path to targets file.
#' @param dmzTargets (numeric; optional) Half-width of m/z interval
#'   centered on target reference m/z value.
#' @param test (logical; optional) if TRUE, shodt run to check the
#'   arguments validity and the compression efficiency
#'
#' @return Returns nothing. It is used for its side effects.
#' @export
#'
#' @examples
compressMS = function(file_in,
                      file_out,
                      ms_type  = 'fticr',
                      mzMin    = 70,
                      mzMax    = 250,
                      dmz      = 0.002,
                      test     = FALSE,
                      compMode = 'grid+targets',
                      targets  = NA,
                      dmzTarget = 0.5
) {

  if(ms_type == 'esquire') {
    cat('\n>>> Nothing to do ! <<<\n')

  } else {

    cat('\n>>> Processing',ms_type,'MS in file:',file_in,'\n')

    #  Get first row
    MS0 = data.table::fread(
      file = file_in ,
      header = FALSE ,
      sep = ',' ,
      stringsAsFactors = FALSE,
      nrows = 1,
      data.table = FALSE
    )
    cat('>>> Size of initial MS in memory:',
        format(object.size(MS0),units='Mb'),'\n')

    mySplit = function(x)
      as.numeric(unlist(strsplit(x, split = ' ')))

    # Build filter on 1rst MS, assuming all MS are on same grid
    ## Get mz
    mz = c()
    msLen = ncol(MS0)-8
    dbl = vapply(
      MS0[1, 9:(msLen + 8)],
      FUN = mySplit,
      FUN.VALUE = numeric(2)
    )
    mz = dbl[1, ]
    if(max(mz) < mzMin | min(mz) > mzMax) {
      message('>>> Error: mzMin/mzMax incompatible with mz range of data')
      stop(call. = FALSE)
    }
    ## Filter mz
    mz0 = mz
    selClosest = 1:length(mz)

    if(grepl('grid',compMode)) {
      # Regular grid
      cat('>>> Grid selection\n')
      mz0 = seq(
        round(max(mzMin,min(mz))),
        round(min(mzMax,max(mz))),
        by = dmz)

      # Pointer to  mz values closest to mz0
      # (timing identical to a vapply implementation)
      selClosest = vector(length=length(mz0))
      for(i in seq_along(mz0))
        selClosest[i] = which.min(abs(mz-mz0[i]))

      # Get rid of doubles
      selClosest = unique(selClosest)
    }

    if(grepl('targets',compMode) &
       !any(is.na(targets))) {
      # Targets-based grid
      cat('>>> Targets selection\n')

      mz0 = mz[selClosest]
      sel = c()
      for(targ in targets) {
        sel = c(sel,
                which(mz0 >= targ - dmzTarget &
                        mz0 <= targ + dmzTarget)
        )
      }
      # Filter and remove doubles (possible targets overlap)
      selClosest = selClosest[unique(sel)]
    }

    # Reorder to compensate for targets +/- random order
    selClosest = sort(selClosest)

    cat('>>> Compress m/z from',length(mz),'to',length(selClosest),
        '(',round(100*(1-length(selClosest)/length(mz))),'%)\n' )

    newRange = range(mz[selClosest])
    cat('>>> New range:',paste0(newRange,collapse='-'),'\n')

    notFinished = TRUE
    j=0
    while(notFinished) {
      j = j+1
      MS1 = c(
        MS0[1,1:6],
        paste0(newRange,collapse ='-'),
        length(selClosest),
        MS0[1,selClosest + 8]
      )
      if(j==1)
        cat('>>> Size of filtered MS in memory:',
            format(object.size(MS1),units='Mb'),'\n')
      if(test) {
        message('>>> Stopped by test=TRUE \n')
        stop(call. = FALSE)
      }
      ## Save to file
      data.table::fwrite(
        MS1,
        file = file_out,
        append = j > 1,
        col.names = FALSE
      )
      cat('+')
      if(j%%20 == 0)
        cat('<',j,'\n')
      MS0 = try(
        # Get new row
        data.table::fread(
          file = file_in ,
          header = FALSE ,
          sep = ',' ,
          stringsAsFactors = FALSE,
          nrows = 1,
          skip  = j,
          data.table = FALSE
        ),
        silent = TRUE
      )
      if(class(MS0) == 'try-error')
        notFinished = FALSE
    }
    cat('\n>>> Finished <<<\n')
  }
}
