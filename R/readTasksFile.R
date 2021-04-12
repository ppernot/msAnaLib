#' Read tasks file
#'
#' @param file (string) a filename
#'
#' @return A dataframe
#'
#' @export
#'
#' @examples
readTasksFile <- function(file) {
  utils::read.table(
    file = file,
    header = TRUE,
    sep = ',',
    check.names = FALSE
  )
}
