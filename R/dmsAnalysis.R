#' Analysis of DMS file
#'
#' @param ms_type
#' @param taskTable
#' @param tgTable
#' @param dataRepo
#' @param figRepo
#' @param tabRepo
#' @param fit_dim
#' @param filter_results
#' @param userTag
#' @param save_figures
#' @param plot_maps
#' @param fallback
#' @param correct_overlap
#' @param weighted_fit
#' @param refine_CV0
#' @param debug
#'
#' @return
#' @export
#'
#' @examples
dmsAnalysis = function(
  ms_type         = c('esquire','fticr')[1],
  taskTable       = NA,
  tgTable         = NA,
  dataRepo        = '../data/',
  figRepo         = '../results/figs/',
  tabRepo         = '../results/tables/',
  fit_dim         = 1,
  filter_results  = TRUE,
  userTag         = paste0('fit_dim_',fit_dim),
  save_figures    = TRUE,
  plot_maps       = FALSE,
  fallback        = TRUE,
  correct_overlap = FALSE,
  weighted_fit    = FALSE,
  refine_CV0      = TRUE,
  debug           = FALSE
) {
  # Changes ----
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
  #
  #---

  options(
    warn = 0,
    stringsAsFactors = FALSE
  )

  # Graphical params for external and local plots
  gPars = setgPars()
  gParsLoc = gPars
  gParsLoc$cex = 1
  gParsLoc$lwd = 1.5

  # Get apparatus-dependent specifications ----
  peakSpecs   = getPeakSpecs(ms_type)

  fwhm_mz_min = peakSpecs$fwhm_mz_min
  fwhm_mz_max = peakSpecs$fwhm_mz_max
  fwhm_mz_nom = peakSpecs$fwhm_mz_nom
  dmz         = peakSpecs$dmz

  fwhm_cv_min = peakSpecs$fwhm_cv_min
  fwhm_cv_max = peakSpecs$fwhm_cv_max
  fwhm_cv_nom = peakSpecs$fwhm_cv_nom
  dCV         = peakSpecs$dCV

  baseline_cor = peakSpecs$baseline_cor
  area_min     = peakSpecs$area_min

  # Gather run params for reproducibility ----
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

  # Check sanity of parameters ----
  assertive::assert_all_are_existing_files(dataRepo)
  assertive::assert_all_are_existing_files(figRepo)
  assertive::assert_all_are_existing_files(tabRepo)

  file = paste0(dataRepo, tgTable)
  assertive::assert_all_are_existing_files(file)

  file = paste0(dataRepo, taskTable)
  assertive::assert_all_are_existing_files(file)

  assertive::assert_is_numeric(fwhm_mz_min)
  if(!assertive::is_positive(fwhm_mz_min))
    stop(paste0('Erreur: fwhm_mz_min =',
                fwhm_mz_min,' should be positive'))

  assertive::assert_is_numeric(fwhm_mz_max)
  if(!assertive::is_positive(fwhm_mz_max))
    stop(paste0('Erreur: fwhm_mz_max =',
                fwhm_mz_max,' should be positive'))

  assertive::assert_is_numeric(fwhm_cv_min)
  if(!assertive::is_positive(fwhm_cv_min))
    stop(paste0('Erreur: fwhm_cv_min =',
                fwhm_cv_min,' should be positive'))

  assertive::assert_is_numeric(fwhm_cv_max)
  if(!assertive::is_positive(fwhm_cv_max))
    stop(paste0('Erreur: fwhm_cv_max =',
                fwhm_cv_max,' should be positive'))

  assertive::assert_is_numeric(area_min)
  if(!assertive::is_positive(area_min))
    stop(paste0('Erreur: area_min =',area_min
                ,' should be positive'))

  assertive::assert_is_numeric(dmz)
  if(!assertive::is_positive(dmz))
    stop(paste0('Erreur: dmz =',dmz,' should be positive'))

  assertive::assert_is_numeric(dCV)
  if(!assertive::is_positive(dCV))
    stop(paste0('Erreur: dCV =',dCV,' should be positive'))

  # Get targets ----
  targets = msAnaLib::readTargetsFile(paste0(dataRepo, tgTable))
  empty = rep(NA,nrow(targets))
  if(!'CV_ref' %in% colnames(targets))
    targets = cbind(targets,CV_ref=empty)

  # Get tasks list ----
  Tasks = msAnaLib::readTasksFile(paste0(dataRepo, taskTable))

  # Check that files exist before proceeding
  if('path' %in% colnames(Tasks)) {
    files = paste0(dataRepo,Tasks[,'path'],Tasks[,'MS_file'])
  } else {
    files = paste0(dataRepo,Tasks[,'MS_file'])
  }
  assertive::assert_all_are_existing_files(files)

  files = paste0(dataRepo,Tasks[,'DMS_file'])
  assertive::assert_all_are_existing_files(files)

  # Loop over tasks ----
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

    ## Get MS ----
    file = paste0(dataRepo, dataPath, msTable)
    lMS  = msAnaLib::getMS(file, ms_type)
    time = lMS$time
    mz   = lMS$mz
    MS   = lMS$MS
    rm(lMS) # Clean-up memory

    # Get CV ----
    file = paste0(dataRepo, CVTable)
    CV0 = utils::read.table(
      file = file,
      header = FALSE,
      sep = '\t',
      stringsAsFactors = FALSE
    )
    CV = rev(CV0[, 4]) # We want increasing CVs

    ## Ensure CV & MS tables conformity
    t0   = Tasks[task,'t0']
    it   = which.min(abs(time - t0))
    selt = which(time >= time[it])
    nt   = length(selt)

    CV0   = Tasks[task,'CV0']
    iCV   = which.min(abs(CV - CV0))
    selCV = which(CV <= CV[iCV])
    nCV   = length(selCV)

    ncut  = min(nt,nCV)
    selt  = selt[1:ncut]
    selCV = rev(rev(selCV)[1:ncut])

    time = time[selt]
    MS   = MS[selt,]
    MS   = apply(MS, 2, rev) # reverse column to conform with CV
    CV   = CV[selCV]
    nCV  = length(CV)

    # Baseline correction ----
    MS = msAnaLib::bslCorMS(MS, baseline_cor)

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

    # Loop over targets ----
    for( it in 1:nrow(targets) ) {

      mz0 = targets[it,'m/z_ref']
      CV0 = targets[it,'CV_ref']

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
          mz0, CV0,
          dmz, dCV,
          mz, CV, MS,
          weighted = weighted_fit,
          fwhm_mz_nom = fwhm_mz_nom
        )
        dimfit = 0
      }

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

        # Quality control
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

      # Plot data and fit results
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
        png(filename = paste0(figRepo, tag, '_', targets[it, 1], '.png'),
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

      # Save XIC and fit file
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

      # if(debug) break
    }

    # Global Heat maps
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

    # Save results
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

    if(debug){
      message('Ended prematurely (debug)...')
      stop(call. = FALSE)
    }

    # End of targets loop ----
  }

  # End of tasks loop ----
}
