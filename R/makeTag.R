#' Create tag
#'
#' @param CVTable
#' @param msTable
#' @param userTag
#'
#' @return
#' @export
#'
#' @examples
makeTag <- function(CVTable, msTable, userTag) {
  date =
    strsplit(
      strsplit(
        CVTable,
        split = ' '
      )[[1]][2],
      split = '-'
    )[[1]][1]

  tag = paste0(
    date,'_',
    strsplit(msTable, split='\\.')[[1]][1],
    ifelse(userTag=='','','_'),userTag
  )
}
