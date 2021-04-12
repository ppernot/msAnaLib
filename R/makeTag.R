#' Create a tag from DMS and CV filenames and a user tag.
#'
#' @param CVTable (string) name of CV file.
#' @param msTable (string) name of ms file.
#' @param userTag (string) string to append to results filenames.
#'
#' @return
#' @export
#'
#' @examples
makeTag <- function(
  CVTable,
  msTable,
  userTag
) {
  # Extract date from CV filename
  date =
    strsplit(
      strsplit(
        CVTable,
        split = ' '
      )[[1]][2],
      split = '-'
    )[[1]][1]
  # Build tag
  tag = paste0(
    date,'_',
    strsplit(msTable, split='\\.')[[1]][1],
    ifelse(userTag=='','','_'),userTag
  )
}
