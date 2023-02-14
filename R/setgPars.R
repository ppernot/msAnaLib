#' Set default graphical parameters
#'
#' @return A list with the following arguments:
#' \describe{
#'   \item{cols}{A set of colors}
#'   \item{cols_tr}{A set of highly transparent colors}
#'   \item{cols_tr2}{Sames as cols_tr with more density}
#'   \item{pty}{A character specifying the type of plot region to be used;
#'     "s" generates a square plotting region and "m" generates the maximal
#'     plotting region (default: "s")}
#'   \item{mar}{Margins around the plot. The default is c(3,3,3,.5).
#'     See \link{graphics::par}.
#'   \item{mgp}{The margin line (in mex units) for the axis title,
#'     axis labels and axis line. The default is c(2,.75,0).}
#'   \item{tcl}{The length of tick marks as a fraction of the height
#'     of a line of text. The default is -0.5.}
#'   \item{lwd}{Width of the plotted lines. The default is 4.}
#'   \item{cex}{A numerical value giving the amount by which
#'     plotting text and symbols should be magnified relative
#'     to the default. The default value is 4.}
#'   \item{cex.leg}{A numerical value giving the amount by which
#'     the legends should be magnified relative
#'     to the default. The default value is 0.7.}
#'   \item{reso}{A numerical value giving the nominal pixel size
#'     for a PNG plot. The default is 1200.}
#' }
#'
#' @export
#'
#' @examples
setgPars = function() {

  myGetColors = function (n,alpha = NULL) {
    color <- c("#E8ECFB","#DDD8EF","#D1C1E1","#C3A8D1","#B58FC2","#A778B4","#9B62A7","#8C4E99","#6F4C9B","#6059A9","#5568B8","#4E79C5","#4D8AC6","#4E96BC","#549EB3","#59A5A9","#60AB9E","#69B190","#77B77D","#8CBC68","#A6BE54","#BEBC48","#D1B541","#DDAA3C","#E49C39","#E78C35","#E67932","#E4632D","#DF4828","#DA2222","#B8221E","#95211B","#721E17","#521A13")
    names(color) <- paste0('color',1:length(color))
    value = seq_along(color)
    value <- scales::rescale(value)
    x <- seq.int(0, 1, length.out = 255)
    color <- (scales::gradient_n_pal(color, values = value))(x)
    pal <- (grDevices::colorRampPalette(color, bias = 1, space = "Lab"))(n)
    if (!is.null(alpha))
      pal <- grDevices::adjustcolor(pal, alpha.f = alpha)
    return(pal)
  }

  list(
    cols     = rev(myGetColors(8))[1:7],
    cols_tr  = rev(myGetColors(8, alpha = 0.2))[1:7],
    cols_tr2 = rev(myGetColors(8, alpha = 0.5))[1:7],
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
