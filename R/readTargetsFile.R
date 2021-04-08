#' Read targets file
#'
#' @param file
#'
#' @return
#' @export
#'
#' @examples
readTargetsFile <- function(file) {
  M = utils::read.table(
    file = file,
    header = TRUE,
    sep = ',',
    dec = '.',
    check.names = FALSE,
    fill = TRUE
  )
  # Backwards separator compatibility
  if(ncol(M) < 4)
    M = utils::read.table(
      file = file,
      header = TRUE,
      sep = ';',
      dec = '.',
      check.names = FALSE,
      fill = TRUE
    )
  # Backwards colnames compatibility
  if('m/z_exact' %in% colnames(M))
    colnames(M)[colnames(M)=='m/z_exact'] = 'm/z_ref'

  return(M)
}
