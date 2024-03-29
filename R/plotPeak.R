#' Plot peak fit results, with two components:
#' (a) localization on the DMS map; (b) 1D peak fit along
#' CV if type == 'CV' or along m/z otherwise.
#'
#' @param mz (vector) m/z vector
#' @param CV (vector) CV vector
#' @param MS (matrix) DMS matrix
#' @param fitOut (object) results of a fit_xD function
#' @param mex (nuleric) exact m/z
#' @param leg (string) legend
#' @param val (string) String describing best-fit parameters
#' @param tag (string) Tag (task+target+user)
#' @param type (string) Which coordinate for right-hand plot (default = 'CV')
#' @param CV0 (numeric) Position of the reference CV value.
#' @param gPars (list) Graphical parameters.
#'
#' @return Nothing. The function is used for its side effects.
#' @export
#'
#' @examples
plotPeak = function(
  mz, CV, MS,
  fitOut,
  mex = NA,
  leg = NA,
  val = NA,
  tag = NA,
  type = 'CV',
  CV0 = NA,
  expandFactor = 5,
  gPars
) {

  nCV = length(CV)

  # Expose fitOut list
  for (n in names(fitOut))
    assign(n,rlist::list.extract(fitOut,n))

  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  par(
    mfrow = c(1, 2),
    mar   = mar,
    mgp   = mgp,
    tcl   = tcl,
    lwd   = lwd,
    cex   = cex
  )

  if(type == 'CV') {
    CVp = CV
    MSp = MS
  }else{
    # Fill matrix CV-wise for better image
    CVp = seq(min(CV),max(CV),length.out=expandFactor*length(CV))
    MSp = matrix(0,nrow=length(CVp),ncol=ncol(MS))
    for(i in 1:length(CV)) {
      iCV = which.min(abs(CVp-CV[i])) # Closest point
      CVp[iCV] = CV[i]
      MSp[iCV,] = MS[i,]
    }
  }

  mzlim  = c(mz1,mz2)
  CVlim  = range(CV)
  if(type == 'CV')
    CVlimf = range(CVf)
  else
    CVlimf = range(CVp)

  ## 1. image
  sel1 = apply(cbind(CVp-CVlim[1],CVp-CVlim[2]),1,prod) <= 0
  sel2 = apply(cbind(mz-mzlim[1],mz-mzlim[2]),1,prod) <= 0

  image(
    CVp[sel1], mz[sel2], MSp[sel1,sel2],
    xlim = CVlim,
    xlab = 'CV',
    ylim = mzlim,
    ylab = 'm/z'
  )
  grid()

  if(class(res) != 'try-error') {
    v   = summary(res)$parameters[,"Estimate"]
    if(length(v) > 3) {
      # Fit 2D
      p = matrix(
        predict(res),
        ncol = ncol(MSloc),
        nrow = nrow(MSloc),
        byrow = TRUE
      )
      contour(CVf, mzloc, p, col = cols_tr2[5], add=TRUE)
    }
    if(type == 'CV')
      rect(CVlimf[1], mz1, CVlimf[2], mz2,
           col = cols_tr[4], border = NA)
  } else {
    v = NA
  }

  if (!is.na(leg))
    title(main = leg, line = 0.5)

  if (!is.na(tag))
    mtext(
      text = tag,
      side = 3,
      line = 2,
      cex = cex
    )

  if (!is.na(CV0))
    abline(v = CV0, lty = 1, col = cols[2])

  if (!is.na(mex))
    abline(h = mex, lty = 1, col = cols[2])

  if (type != 'CV' & !any(is.na(v))) {
    abline(h = v[1], lty = 2, col = cols[2])
  } else {
    if( length(v) == 3) {
      abline(v = v[1], lty = 2, col = cols[2])
    } else {
      abline(h = v[1], lty = 2, col = cols[2])
      abline(v = v[3], lty = 2, col = cols[2])
    }
    if(!is.na(mz0))
      abline(h = mz0, lty = 2, col = cols[2])
  }

  if(type == 'CV') {
    ## 2. CV profile
    xmod = seq(min(CV),max(CV),length.out = 1000)
    if(class(res) != 'try-error') {
      v   = summary(res)$parameters[,"Estimate"]
      vmod = peakShape(xmod, v)
    } else {
      v = NA
      vmod = rep(NA,length(xmod))
    }
    ylim = range(c(mMS,vmod), na.rm = TRUE, finite = TRUE)
    plot(
      CV, mMStot,
      type = 'l', #pch = 16,
      col  = cols[4],
      xlim = CVlim,
      xlab = 'CV',
      ylim = ylim,
      ylab = 'a.u.'
    )
    title(main = 'Integ. CV profile', line = 0.5)
    rect(CVlimf[1], -100, CVlimf[2], ylim[2]*1.2,
         col = cols_tr[4], border=NA)

    ## 2.1 Gaussian fit
    if(!any(is.na(v))) {
      lines(xmod, vmod, col = cols[6])
      if(length(v) == 3)
        abline(v = v[1], lty = 2, col = cols[2])
      else
        abline(v = v[3], lty = 2, col = cols[2])
    }

    if (!is.na(CV0))
      abline(v = CV0, lty = 1, col = cols[2])

    ## Add fit results
    if (!is.na(val))
      legend(
        x=CVlim[1],
        y=ylim[2]-0.75*strheight(val,units='user',cex=0.8),
        yjust = 0,
        # title  = val,
        legend = val,
        bty    = 'n',
        cex    = 0.8)

  } else {

    ## 2. MS profile
    xmod = seq(mz1,mz2,length.out = 1000)
    if(class(res) != 'try-error') {
      v   = summary(res)$parameters[,"Estimate"]
      vmod = peakShape(xmod, v)
    } else {
      v = NA
      vmod = rep(NA,length(xmod))
    }
    ylim = c(
      0,
      1.2*max(c(MSl,vmod), na.rm = TRUE, finite = TRUE)
    )
    plot(
      mzl, MSl,
      type = 'l', #pch = 16,
      col  = cols[4],
      xlim = mzlim,
      xlab = 'm/z',
      ylim = ylim,
      ylab = 'a.u.',
      yaxs = 'i'
    )
    title(main = 'MS profile', line = 0.5)

    ## 2.1 Gaussian fit
    if(!any(is.na(v))) {
      lines(xmod, vmod, col = cols[6])
      abline(v = v[1], col = cols[2], lty = 2)
    }

    if (!is.na(mex))
      abline(v = mex, lty = 1, col = cols[2])

    ## Add fit results
    if (!is.na(val))
      legend(
        x = mzlim[1],
        y = ylim[2] - strheight(val, units = 'user', cex = 0.8),
        yjust = 0,
        legend = val,
        bty    = 'n',
        cex    = 0.8)
  }

}
