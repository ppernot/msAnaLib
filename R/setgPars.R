#' Set graphical parameters
#'
#' @return
#' @export
#'
#' @examples
setgPars = function() {
  list(
    cols     = rev(inlmisc::GetColors(8))[1:7],
    cols_tr  = rev(inlmisc::GetColors(8, alpha = 0.2))[1:7],
    cols_tr2 = rev(inlmisc::GetColors(8, alpha = 0.5))[1:7],
    pty      = 's',
    mar      = c(3,3,3,.5),
    mgp      = c(2,.75,0),
    tcl      = -0.5,
    lwd      = 4.0,
    cex      = 4.0,
    cex.leg  = 0.7,
    reso     = 1200
  )
}
