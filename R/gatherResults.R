#' Gather results in data frame
#'
#' @param Tasks
#' @param tabRepo
#' @param userTag
#'
#' @return
#' @export
#'
#' @examples
gatherResults <- function(Tasks, tabRepo, userTag) {
  D = NULL
  for(task in 1:nrow(Tasks)) {

    msTable = Tasks[task,'MS_file']
    CVTable = Tasks[task,'DMS_file']
    # dilu    = Tasks[task,'dilu']

    # Build tag
    tag = makeTag(CVTable, msTable, userTag)
    file = paste0(tabRepo,tag,'_results.csv')
    if(!file.exists(file))
      stop(paste0('Missing file:',file))

    M = read.csv(file = file, check.names = FALSE)
    # M = cbind(M,dilu)
    D = rbind(D,M)
  }
  return(D)
}
