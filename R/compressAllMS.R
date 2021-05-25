#' Compress a set of DMS files.
#'
#' @param dataRepo (string; optional) Path to the data files.
#' @param origMsDir (string; mandatory) Path to repertory of DMS files
#'   to be compressed
#' @param compMsDir (string; mandatory) Path to repertory of compressed files
#' @param ms_type (string; optional) Type of mass spectrum,
#'   one of `esquire` or `fticr` (default).
#'   Compression only for `fticr` MS.
#' @param compMode (string; optional) Compression mode:
#'   either `grid`, `targets` (default),
#'   or their combination `grid+targets`.
#' @param mzMin (numeric; optional) Min m/z value of compression grid
#' @param mzMax (numeric; optional) Max m/z value of compression grid
#' @param dmz (numeric; optional) m/z step or compression grid
#' @param tgTable (string; optional)  Path to targets file.
#' @param dmzTarget (numeric; optional) Half-width of m/z interval
#'   centered on target reference m/z value.
#' @param test (logical; optional) if TRUE, shodt run to check the
#'   arguments validity and the compression efficiency
#'
#' @return Returns nothing. It is used for its side effects.
#' @export
#'
#' @examples
compressAllMS <- function(
  origMsDir ,
  compMsDir ,
  ms_type   = c('fticr','esquire'),
  dataRepo  = '../data/',
  compMode  = c('targets','grid','grid+targets'),
  tgTable   = NA,
  dmzTarget = 0.5,
  mzMin     = 70,
  mzMax     = 250,
  dmz       = 0.001,
  test      = FALSE
) {

  options(stringsAsFactors = FALSE)

  ms_type  = match.arg(ms_type)
  compMode = match.arg(compMode)

  # Checks and Get ms list ####
  assertive::assert_all_are_existing_files(dataRepo)

  origMs = file.path(dataRepo, origMsDir)
  assertive::assert_all_are_existing_files(origMs)

  msFiles = list.files(origMs) # All files in source Dir
  if (length(msFiles) == 0) {
    message('>>> No data files to process !!!')
    stop(call. = FALSE)
  }

  compMs = file.path(dataRepo, compMsDir)
  if (!dir.exists(compMs))
    dir.create(compMs)

  targets = NA
  if(grepl('targets',compMode)) {
    file = file.path(dataRepo, tgTable)
    assertive::assert_all_are_existing_files(file)
    targets = readTargetsFile(file)[,'m/z_ref']

    assertive::assert_is_numeric(dmzTarget)
    if(!assertive::is_positive(dmzTarget)) {
      message('>>> Error: dmzTarget =', dmzTarget,' should be positive')
      stop(call. = FALSE)
    }
  }

  if(grepl('grid',compMode)) {
    for(arg in c('mzMin','mzMax','dmz')) {
      val = get(arg)
      assertive::assert_is_numeric(val)
      if(!assertive::is_positive(val)) {
        message('>>> Error: ',arg,' =', val,' should be positive')
        stop(call. = FALSE)
      }
    }
  }

  # Loop over files ####
  for (ms in msFiles)
    compressMS(
      file_in   = file.path(origMs, ms),
      file_out  = file.path(compMs, ms),
      ms_type   = ms_type,
      mzMin     = mzMin,
      mzMax     = mzMax,
      dmz       = dmz,
      test      = test,
      compMode  = compMode,
      targets   = targets,
      dmzTarget = dmzTarget
    )
}
