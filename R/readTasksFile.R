#' Read tasks file
#'
#' @param file
#'
#' @return
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
