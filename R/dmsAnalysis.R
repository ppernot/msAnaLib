#' Analysis of DMS file
#'
#' @param ms_type (string; optional) Type of mass spectrum,
#'   one of `esquire` (default) or `fticr`.
#' @param taskTable (string; mandatoty) Path to tasks file.
#' @param tgTable (string; mandatoty) Path to targets file.
#' @param dataRepo (string; optional) Path to the data files.
#' @param figRepo (string; optional) Path to output figures.
#' @param tabRepo (string; optional) Path to output tables.
#' @param fit_dim (numeric; optional) Peak fit options.
#'   0: 1D fit along m/z for data at given CV value;
#'   1 (default): 1D fit along CV for data averaged in m/z interval;
#'   2: 2D fit.
#' @param filter_results (logical, optional) Flag to filter
#'   fit results according to their fwhm and area.
#'   Constraints are provided in `peakSpecs` list.
#' @param userTag (string; optional) string to append to results
#'   filenames.
#' @param save_figures (logical; optional) save figures to png files
#'   (default: TRUE).
#' @param plot_maps (logical; optional) plot and save DMS maps
#'   to png files (default: FALSE).
#' @param fallback (logical; optional) Fallback on 1D fit
#'   (\code{fit_dim = 1}) when 2D fit (\code{fit_dim = 2})
#'   fails (default: TRUE).
#' @param correct_overlap (logical; optional) Use a correction
#'   for the peaks overlap (logical: FALSE). Experimental, to be
#'   implemented.
#' @param weighted_fit (logical; optional) Used weighted data for
#'   the fit (default: FALSE).
#' @param refine_CV0 (logical; optional) Refine the CV window to
#'   improve peak fit (default: TRUE).
#' @param debug (logical; optional) Stop after the fitst task
#'   (default: FALSE).
#' @param gPars (list; optional) Supersedes the default graphical
#'   parameters values. See \link{setgPars}.
#' @param peakSpecs (list; optional) Supersedes the default peak
#'   specifications. See \link{getPeakSpecs}.
#'   values.
#'
#' @return The function is mostly used for its side effects.
#'   It invisibly returns a dataframe of results issued from
#'   \link{gatherResults}.
#'
#' @export
#'
#' @examples
dmsAnalysis = function(
  taskTable,
  tgTable,
  dataRepo        = '../data/',
  figRepo         = '../results/figs/',
  tabRepo         = '../results/tables/',
  ms_type         = c('esquire','fticr'),
  fit_dim         = c(1,2,0),
  filter_results  = TRUE,
  userTag         = paste0('fit_dim_',fit_dim),
  save_figures    = TRUE,
  plot_maps       = FALSE,
  fallback        = TRUE,
  correct_overlap = FALSE,
  weighted_fit    = FALSE,
  refine_CV0      = TRUE,
  debug           = FALSE,
  gPars           = list(),
  peakSpecs       = list()
) {
  #*********************************************************
  # Changes ----
  #*********************************************************
  #
  # 2020_07_16 [PP]
  # - Added 2 new colums to results table
  #   to store m/z and u_m/z from 2D fit
  # - Integrate fast method (fit_dim = 0)
  # - Added data path management
  # - Added optional tag ('userTag'):
  #   presently defined by 'fit_dim' value,
  #   to avoid overwriting of results files
  #   when trying different 'fit_dim' options
  # 2020_07_20 [PP]
  # - Replaced '=' by '_' in userTag (Windows pb.)
  # - Reparameterized gaussians with area replacing height
  #   (avoids covariances in estimation of u_area )
  # - Suppressed 'rho' param in 2D gaussians
  # 2020_07_21 [PP]
  # - Corrected typo in calculation of u_area in 'getPars1D()'
  # 2020_07_24 [PP]
  # - Save ctrl params as metadata
  # 2020_09_24 [PP]
  # - All input files should now be comma-delimited
  # 2021_03_15 [PP]
  # - Adapted code to FT-ICR MS files:
  #   * new 'ms_type' flag
  #   * new 'getMS()' function to handle ms_type cases
  #   * new 'trapz()' function to handle MS integration
  #     with irregular mz grids
  #   * adapted 'fit1D()' and 'fit2D()' to use 'trapz()'
  # 2021_04_06 [PP]
  # - Added baseline correction function
  #   * new 'baseline_cor' flag
  #   * new 'bslCorMS()' function
  # 2021_04_07 [PP]
  # - Changed management of peak specifications to facilitate
  #   modifications of apparatus characteristics
  #   * new 'getPeakSpecs()' functions
  #   * changed 'const_fwhm' to 'fwhm_cv_nom'
  #   * new 'fwhm_mz_nom' variable
  #   * changed logic on fwhm fit constraints
  #     Fit is now always constrained in fit1D_MS(),
  #     fit_1D() and fit_2D()
  # 2021_04_09 [PP]
  # - Changed to function in msAnaLib package
  #
  #*********************************************************

  options(stringsAsFactors = FALSE)

  #*********************************************************
  # Graphical params for external and local plots
  #*********************************************************
  ## Get default params
  gParsDef = setgPars()
  ## Override by user's specs, if any
  if(!is.list(gPars)) {
    message('>>> Argument gPars should be a list !\n')
    stop(call.=FALSE)
  } else {
    if(length(gPars) != 0) {
      for (n in names(gPars)) {
        if(!n %in% names(gParsDef)) {
          message('>>> Incorrect variable name in gPars: ',n,'\n')
          stop(call.=FALSE)
        }
        gParsDef[[n]] = rlist::list.extract(gPars, n)
      }
    }
  }
  gPars =gParsDef
  gParsLoc = gPars
  gParsLoc$cex = 1
  gParsLoc$lwd = 1.5

  #*********************************************************
  # Set apparatus-dependent specifications ----
  #*********************************************************
  ## Get default specs
  peakSpecsDef = getPeakSpecs(ms_type)
  for (n in names(peakSpecsDef))
    assign(n, rlist::list.extract(peakSpecsDef, n))
  ## Override by user's specs, if any
  if(!is.list(peakSpecs)) {
    message('>>> Argument peakSpecs should be a list !\n')
    stop(call.=FALSE)
  } else {
    if(length(peakSpecs) != 0) {
      for (n in names(peakSpecs)) {
        if(!n %in% names(peakSpecsDef)) {
          message('>>> Incorrect variable name in peakSpecs: ',n,'\n')
          stop(call.=FALSE)
        }
        assign(n, rlist::list.extract(peakSpecs, n))
      }
    }
  }
  #*********************************************************

  #*********************************************************
  # Check sanity of parameters ----
  #*********************************************************

  ms_type = match.arg(ms_type)
  if(!fit_dim %in% 0:2) {
    message('>>> Error: fit_dim should be 0, 1, or 2 !')
    stop(call. = FALSE)
  }

  assertive::assert_all_are_existing_files(dataRepo)
  assertive::assert_all_are_existing_files(figRepo)
  assertive::assert_all_are_existing_files(tabRepo)

  file = file.path(dataRepo, tgTable)
  assertive::assert_all_are_existing_files(file)

  file = file.path(dataRepo, taskTable)
  assertive::assert_all_are_existing_files(file)

  for(arg in c('fwhm_mz_min','fwhm_mz_max',
               'fwhm_cv_min','fwhm_cv_max',
               'area_min','dmz','dCV')     ) {
    val = get(arg)
    assertive::assert_is_numeric(val)
    if(!assertive::is_positive(val)) {
      message('>>> Error: ',arg,' =', val,' should be positive')
      stop(call. = FALSE)
    }
  }
  #*********************************************************

  #*********************************************************
  # Gather run params for reproducibility ----
  #*********************************************************
  ctrlParams = list(
    userTag        = userTag,
    ms_type        = ms_type,
    taskTable      = taskTable,
    tgTable        = tgTable,
    dataRepo       = dataRepo,
    filter_results = filter_results,
    fwhm_mz_min    = fwhm_mz_min,
    fwhm_mz_max    = fwhm_mz_max,
    fwhm_mz_nom    = fwhm_mz_nom,
    dmz            = dmz,
    fwhm_cv_min    = fwhm_cv_min,
    fwhm_cv_max    = fwhm_cv_max,
    fwhm_cv_nom    = fwhm_cv_nom,
    dCV            = dCV,
    area_min       = area_min,
    fit_dim        = fit_dim,
    fallback       = fallback,
    weighted_fit   = weighted_fit,
    refine_CV0     = refine_CV0,
    baseline_cor   = baseline_cor
  )
  #*********************************************************

  #*********************************************************
  # Get targets ----
  #*********************************************************
  targets = msAnaLib::readTargetsFile(paste0(dataRepo, tgTable))
  empty = rep(NA, nrow(targets))
  if (!'CV_ref' %in% colnames(targets))
    targets = cbind(targets, CV_ref = empty)
  #*********************************************************

  #*********************************************************
  # Get tasks list ----
  #*********************************************************
  Tasks = msAnaLib::readTasksFile(paste0(dataRepo, taskTable))
  #*********************************************************

  #*********************************************************
  # Check that files exist before proceeding
  #*********************************************************
  if('path' %in% colnames(Tasks)) {
    files = paste0(dataRepo, Tasks[, 'path'], Tasks[, 'MS_file'])
  } else {
    files = paste0(dataRepo, Tasks[, 'MS_file'])
  }
  assertive::assert_all_are_existing_files(files)

  files = paste0(dataRepo,Tasks[,'DMS_file'])
  assertive::assert_all_are_existing_files(files)
  #*********************************************************

  #*********************************************************
  # Loop over tasks ----
  #*********************************************************
  dilu = NA
  for(task in 1:nrow(Tasks)) {

    msTable = Tasks[task,'MS_file']
    CVTable = Tasks[task,'DMS_file']
    if('dilu' %in% colnames(Tasks))
      dilu    = Tasks[task,'dilu']
    dataPath = ''
    if('path' %in% colnames(Tasks))
      if(!is.na(Tasks[task,'path']))
        dataPath = Tasks[task,'path']

    # Build tag
    tag  = msAnaLib::makeTag(CVTable, msTable, userTag)

    #*********************************************************
    ## Get MS ----
    #*********************************************************
    file = paste0(dataRepo, dataPath, msTable)
    lMS  = msAnaLib::getMS(file, ms_type)
    time = lMS$time
    mz   = lMS$mz
    MS   = lMS$MS
    rm(lMS) # Clean-up memory
    #*********************************************************

    #*********************************************************
    # Get CV ----
    #*********************************************************
    file = paste0(dataRepo, CVTable)
    CV0 = utils::read.table(
      file = file,
      header = FALSE,
      sep = '\t',
      stringsAsFactors = FALSE
    )
    CV = CV0[, 4]
    # We want increasing CVs
    reverseCV = FALSE
    if(CV[2] < CV[1])
      reverseCV = TRUE
    if(reverseCV)
      CV = rev(CV)
    #*********************************************************

    #*********************************************************
    ## Ensure CV & MS tables conformity
    #*********************************************************
    t0   = Tasks[task,'t0']
    it   = which.min(abs(time - t0))
    selt = which(time >= time[it])
    nt   = length(selt)

    CV0   = Tasks[task,'CV0']
    iCV   = which.min(abs(CV - CV0))
    selCV = which(CV >= CV[iCV])
    if(reverseCV)
      selCV = which(CV <= CV[iCV])
    nCV   = length(selCV)

    ncut  = min(nt,nCV)
    selt  = selt[1:ncut]
    selCV = selCV[1:ncut]
    if(reverseCV)
      selCV = rev(rev(selCV)[1:ncut])

    time = time[selt]
    MS   = MS[selt,]
    if(reverseCV)
      MS   = apply(MS, 2, rev) # reverse column to conform with CV
    CV   = CV[selCV]
    nCV  = length(CV)
    #*********************************************************

    #*********************************************************
    # Baseline correction ----
    #*********************************************************
    MS = msAnaLib::bslCorMS(MS, baseline_cor)
    #*********************************************************

    ## Initialize results table
    resu = cbind(
      targets,
      empty,empty,
      empty,empty,
      empty,empty,
      empty,empty,
      empty,empty,
      empty,empty,
      empty)
    colnames(resu) = c(
      colnames(targets),
      'm/z',     'u_m/z',
      'CV',      'u_CV',
      'FWHM_m/z','u_FWHM_m/z',
      'FWHM_CV', 'u_FWHM_CV',
      'Area',    'u_Area',
      'fit_dim',  'dilu',
      'tag'
    )

    if( fit_dim == 0) {
      xic = matrix(mz,ncol=1)
      colnames(xic) = 'm/z'
      xfi = matrix(mz,ncol=1)
      colnames(xfi) = 'm/z'
    } else {
      xic = cbind(time,rev(CV))
      colnames(xic) = c('time','CV')
      xfi = cbind(time,rev(CV))
      colnames(xfi) = c('time','CV')
    }

    #*********************************************************
    # Loop over targets ----
    #*********************************************************
    for( it in 1:nrow(targets) ) {

      mz0 = targets[it,'m/z_ref']
      CV0 = targets[it,'CV_ref']

      #*********************************************************
      ## Fit ----
      #*********************************************************
      if(fit_dim == 2) {
        # 2D fit of peaks
        fitOut = msAnaLib::fit2D(
          mz0, CV0,
          dmz, dCV,
          mz, CV, MS,
          fwhm_mz_nom = fwhm_mz_nom,
          fwhm_cv_nom = fwhm_cv_nom,
          weighted = weighted_fit,
          refine_CV0 = refine_CV0,
          correct_overlap = correct_overlap
        )
        dimfit = 2
        if(class(fitOut$res) == 'try-error' & fallback) {
          # 1D fit of peaks
          fitOut = msAnaLib::fit1D(
            mz0, CV0,
            dmz, dCV,
            mz, CV, MS,
            fwhm_cv_nom = fwhm_cv_nom,
            weighted = weighted_fit,
            refine_CV0 = refine_CV0,
            correct_overlap = correct_overlap
          )
          dimfit = 1
        }

      } else if (fit_dim == 1) {
        # 1D fit of peaks; fixed m/z
        fitOut = msAnaLib::fit1D(
          mz0, CV0,
          dmz, dCV,
          mz, CV, MS,
          fwhm_cv_nom = fwhm_cv_nom,
          weighted = weighted_fit,
          refine_CV0 = refine_CV0,
          correct_overlap = correct_overlap
        )
        dimfit = 1

      } else {
        # 1D m/z fit of peak; fixed CV
        fitOut = msAnaLib::fit1D_MS(
          mz0, CV0, dmz,
          mz, CV, MS,
          weighted = weighted_fit,
          fwhm_mz_nom = fwhm_mz_nom
        )
        dimfit = 0
      }
      #*********************************************************

      #*********************************************************
      # Extract fit results ----
      #*********************************************************
      for (n in names(fitOut))
        assign(n,rlist::list.extract(fitOut,n))

      targets[it,'m/z_ref'] = mz0

      if(class(res)=="try-error") {
        # Fit failed => no fit params
        v       = NA
        mzopt   = NA
        cvopt   = NA
        fwhm_mz = NA
        fwhm_cv = NA
        area    = NA
        warning = TRUE

      } else {
        v   = summary(res)$parameters[,"Estimate"]
        peakPars = getPars(res,dimfit)
        for (n in names(peakPars))
          assign(n,rlist::list.extract(peakPars,n))

        #*********************************************************
        # Quality control
        #*********************************************************
        if(filter_results &
           (
             ifelse(
               !is.na(fwhm_cv),
               fwhm_cv <= fwhm_cv_min | fwhm_cv >= fwhm_cv_max,
               FALSE
             ) |
             ifelse(
               !is.na(fwhm_mz),
               fwhm_mz <= fwhm_mz_min | fwhm_mz >= fwhm_mz_max,
               FALSE
             ) |
             area <= area_min
           )
        ) {
          warning = TRUE
          # Do not store results

        } else {
          warning = FALSE
          # Store in results table
          resu[it,'m/z']        = signif(mzopt,6)
          resu[it,'u_m/z']      = signif(u_mz,2)
          resu[it,'CV']         = signif(cvopt,4)
          resu[it,'u_CV']       = signif(u_cv,2)
          resu[it,'FWHM_m/z']   = signif(fwhm_mz,3)
          resu[it,'u_FWHM_m/z'] = signif(u_fwhm_mz,2)
          resu[it,'FWHM_CV']    = signif(fwhm_cv,3)
          resu[it,'u_FWHM_CV']  = signif(u_fwhm_cv,2)
          resu[it,'Area']       = signif(area,3)
          resu[it,'u_Area']     = signif(u_area,2)
        }
      }
      resu[it,'fit_dim'] = dimfit
      resu[it,'dilu'] = dilu
      resu[it,'tag'] = tag
      #*********************************************************

      #*********************************************************
      # Plot data and fit results
      #*********************************************************
      pars = paste0(
        ifelse (warning, '** WARNING **\n','') ,
        ifelse(dimfit == 1,
               '',
               paste0('m/z = ', signif(mzopt,6),'\n')),
        ifelse(dimfit == 0,
               '',
               paste0('CV = ', signif(cvopt,4),'\n')),
        ifelse(dimfit ==1,
               '',
               paste0('FWHM_m/z = ', signif(fwhm_mz,3),'\n')),
        ifelse(dimfit ==0,
               '',
               paste0('FWHM_CV = ', signif(fwhm_cv,3),'\n')),
        'Area = ', signif(area,3)
      )

      if (class(fitOut$res) != "try-error")
        print(coefficients(fitOut$res))

      msAnaLib::plotPeak(
        mz, CV, MS,
        fitOut,
        mex = targets[it,'m/z_ref'],
        leg = targets[it,'Name'],
        tag = tag,
        val = pars,
        type = ifelse(dimfit==0,'m/z','CV'),
        CV0 = CV0,
        gPars = gParsLoc
      )

      if(save_figures) {
        png(filename = paste0(figRepo, tag, '_',
                              targets[it, 1], '.png'),
            width    = 2*gPars$reso,
            height   =   gPars$reso )
        msAnaLib::plotPeak(
          mz, CV, MS,
          fitOut,
          mex = targets[it,'m/z_ref'],
          leg = targets[it,'Name'],
          tag = tag,
          val = pars,
          type = ifelse(dimfit==0,'m/z','CV'),
          CV0 = CV0,
          gPars = gPars
        )
        dev.off()
      }
      #*********************************************************

      #*********************************************************
      # Save XIC and fit file
      #*********************************************************
      nam0 = colnames(xic)
      if(fit_dim == 0) {
        xic = cbind(xic,mMStot)
        fit = msAnaLib::peakShape(mz, v)
        xfi = cbind(xfi,fit)
      } else {
        xic = cbind(xic,rev(mMStot))
        fit = msAnaLib::peakShape(CV, v)
        xfi = cbind(xfi,rev(fit))
      }
      colnames(xic) = c(nam0,targets[it,1])
      colnames(xfi) = c(nam0,targets[it,1])
      #*********************************************************

      # Stop after first target in first task...
      if(debug)
        break()
    }

    #*********************************************************
    # End of targets loop ----
    #*********************************************************

    #*********************************************************
    # Global Heat maps
    #*********************************************************
    if(plot_maps) {
      mex = targets[,'m/z_ref']
      msAnaLib::plotMaps(
        mz, CV, MS,
        mex = mex,
        leg = 'log10(MS)',
        tag = tag,
        mzlim = c(min(mex)-5*dmz,max(mex)+5*dmz),
        CVlim = range(CV),
        logz = TRUE,
        gPars = gParsLoc
      )
      if(save_figures) {
        png(
          filename = paste0(figRepo, tag,
                            '_heatmaps.png'),
          width    = 2*gPars$reso,
          height   =   gPars$reso )
        msAnaLib::plotMaps(
          mz, CV, MS,
          mex = mex,
          leg = 'log10(MS)',
          tag = tag,
          mzlim = c(min(mex)-5*dmz,max(mex)+5*dmz),
          CVlim = range(CV),
          logz = TRUE,
          gPars = gPars
        )
        dev.off()
      }
    }
    #*********************************************************

    #*********************************************************
    ## Save results ----
    #*********************************************************
    write.csv(resu,
              file = paste0(tabRepo, tag, '_results.csv'),
              row.names = FALSE)
    write.csv(xic,
              file  = paste0(tabRepo, tag, '_XIC.csv'),
              row.names = FALSE)
    write.csv(xfi,
              file  = paste0(tabRepo, tag, '_fit.csv'),
              row.names = FALSE)
    # Metadata
    rlist::list.save(
      ctrlParams,
      file = file.path(
        tabRepo,
        paste0(tag,'_ctrlParams.yaml')
      )
    )
    #*********************************************************

    if(debug) {
      message('>>> Stopped by debug = TRUE...\n')
      break()
    }

  }

  #*********************************************************
  # End of tasks loop ----
  #*********************************************************

  invisible(resu)
}
