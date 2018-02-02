Function banyan_sigma, stars_data, COLUMN_NAMES=column_names, HYPOTHESES=hypotheses, LN_PRIORS=ln_priors, $
  UNIT_PRIORS=unit_priors, NTARGETS_MAX=ntargets_max, LNP_ONLY=lnp_only, NO_XYZ=no_xyz, $
  RA=ra,DEC=dec,PMRA=pmra,PMDEC=pmdec,EPMRA=epmra,EPMDEC=epmdec,DIST=dist,EDIST=edist,RV=rv,ERV=erv, $
  PSIRA=psira,PSIDEC=psidec,EPSIRA=epsira,EPSIDEC=epsidec, PLX=plx, EPLX=eplx, $
  CONSTRAINT_DIST_PER_HYP=constraint_dist_per_hyp, CONSTRAINT_EDIST_PER_HYP=constraint_edist_per_hyp, $
  USE_RV=use_rv, USE_DIST=use_dist, USE_PLX=use_plx, USE_PSI=use_psi, OVERRIDE_ERRORS=override_errors
  
  ;+
  ; NAME:
  ;       BANYAN_SIGMA
  ;
  ; PURPOSE:
  ;       Calculate the membership probability that a given astrophysical object belongs to one of the currently
  ;       known 27 young associations within 150 pc of the Sun, using Bayesian inference. This tool uses the sky
  ;       position and proper motion measurements of an object, with optional radial velocity (RV) and distance (D)
  ;       measurements, to derive a Bayesian membership probability. By default, the priors are adjusted such that
  ;       a probability treshold of 90% will recover 50%, 68%, 82% or 90% of true association members depending on
  ;       what observables are input (only sky position and proper motion, with RV, with D, with both RV and D,
  ;       respectively).
  ;       
  ;       Please see Gagné et al. 2017 (ApJ, XX, XX) for more detail.
  ;       
  ;       An online version of this tool is available for 1-object queries at http://www.exoplanetes.umontreal.ca/banyan/banyansigma.php.
  ;       
  ; REQUIREMENTS:
  ;       (1) A fits file containing the parameters of the multivariate Gaussian models of each Bayesian hypothesis must be included
  ;           at /data/banyan_sigma_parameters.fits in the directory where BANYAN_SIGMA.pro is compiled. 
  ;           The fits file can be written with MWRFITS and must contain an IDL array of structures of N elements, where N is the total
  ;           number of multivariate Gaussians used in the models of all Bayesian hypotheses. Each element of this structure contains
  ;           the following information:
  ;           NAME - The name of the model (scalar string).
  ;           CENTER_VEC - Central XYZUVW position of the model (6D vector, in units of pc and km/s).
  ;           COVARIANCE_MATRIX - Covariance matrix in XYZUVW associated with the model (6x6 matrix, in mixed units of pc and km/s).
  ;           PRECISION_MATRIX - (Optional) Matrix inverse of COVARIANCE_MATRIX, to avoid re-calculating it many times (6x6 matrix).
  ;           LN_NOBJ - (Optional) Natural logarithm of the number of objects used to build the synthetic model (scalar). This is not used in BANYAN_SIGMA.
  ;           COVARIANCE_DETERM - (Optional) Determinant of the covariance matrix, to avoid re-calculating it many times (scalar).
  ;           PRECISION_DETERM - (Optional) Determinant of the precision matrix, to avoid re-calculating it many times (scalar).
  ;           LN_ALPHA_K - (Optional) Natural logarithm of the alpha_k inflation factors that ensured a fixed rate of true positives at a given Bayesian
  ;                        probability treshold. See Gagné et al. 2017 (ApJ, XX, XX) for more detail (scalar or 4-elements vector).
  ;                        This is not used in BANYAN_SIGMA.
  ;           LN_PRIOR - Natural logarithm of the Bayesian prior (scalar of 4-elements vector). When this is a 4-elements vector, the cases with
  ;                      only proper motion, proper motion + radial velocity, proper motion + distance or proper motion + radial velocity + distance
  ;                      will be used with the corresponding element of the LN_PRIOR vector.
  ;           LN_PRIOR_OBSERVABLES - Scalar string or 4-elements vector describing the observing modes used for each element of LN_PRIOR.
  ;                                  This is not used in BANYAN_SIGMA.
  ;           COEFFICIENT - Coefficient (or weight) for multivariate Gaussian mixture models. This will only be used if more than one element of the
  ;                         parameters array have the same model name (see below).  
  ;           
  ;           When more than one elements have the same model name, BANYAN_SIGMA will use the COEFFICIENTs to merge its Bayesian probability,
  ;           therefore representing the hypothesis with a multivariate Gaussian model mixture.
  ;           
  ;       (2) (Optional) A fits file containing the various performance metrics (true positive rate, false positive rate, positive
  ;           predictive value) as a function of the Bayesian probability treshold, for each young association. Each element of this structure contains
  ;           the following information:
  ;           NAME - The name of the model (scalar string).
  ;           PROBS - N-elements array containing a list of Bayesian probabilities (%).
  ;           TPR - Nx4-elements array containing the rate of true positives that correspond to each of the Bayesian probability (lower) tresholds stored in PROBS.
  ;           FPR - Nx4-elements array containing the rate  of false positives that correspond to each of the Bayesian probability (lower) tresholds stored in PROBS.
  ;           PPV - Nx4-elements array containing the Positive Predictive Values that correspond to each of the Bayesian probability (lower) tresholds stored in PROBS.
  ;           NFP - Number of expected false positives (FPR times the ~7 million stars in the Besancon simulation of the Solar neighborhood) 
  ;           
  ;           Each component of the 4-elements dimension of TPR, FPR, NFP and PPV corresponds to a different mode of input data,
  ;           see the description of "LN_PRIOR" above for more detail.
  ;           
  ;           When this fits file is used, the Bayesian probabilities of each star will be associated with a TPR, FPR, NFP and PPV values in the METRICS sub-structure of
  ;           the output structure.
  ;           
  ;           This file must be located at /data/banyan_sigma_metrics.fits in the directory where BANYAN_SIGMA.pro is compiled.
  ;           
  ; CALLING SEQUENCE:
  ; 
  ;       OUTPUT_STRUCTURE = BANYAN_SIGMA([ stars_data, COLUMN_NAMES=column_names, HYPOTHESES=HYPOTHESES, LN_PRIORS=LN_PRIORS,
  ;           NTARGETS_MAX=ntargets_max, RA=RA, DEC=DEC, PMRA=PMRA, PMDEC=PMDEC, EPMRA=EPMRA, EPMDEC=EPMDEC, DIST=DIST,
  ;           EDIST=EDIST, RV=RV, ERV=ERV, PSIRA=PSIRA, PSIDEC=PSIDEC, EPSIRA=EPSIRA, EPSIDEC=EPSIDEC, PLX=PLX, EPLX=EPLX,
  ;           CONSTRAINT_DIST_PER_HYP=CONSTRAINT_DIST_PER_HYP, CONSTRAINT_EDIST_PER_HYP=CONSTRAINT_EDIST_PER_HYP,
  ;          /UNIT_PRIORS, /LNP_ONLY, /NO_XYZ, /USE_RV, /USE_DIST, /USE_PLX, /USE_PSI ])
  ;
  ; OPTIONAL INPUTS:
  ;       stars_data - An IDL structure (or array of structures when more than one objects are analyzed) that contain at least the following tags:
  ;                    RA, DEC, PMRA, PMDEC, EPMRA, and EPMDEC. It can also optionally contain the tags RV, ERV, DIST, EDIST, PLX, EPLX
  ;                    PSIRA, PSIDEC, EPSIRA, EPSIDEC. See the corresponding keyword descriptions for more information.
  ;                    If this input is not used, the keywords RA, DEC, PMRA, PMDEC, EPMRA, and EPMDEC must all be specified.
  ;       column_names - An IDL structure that contains the names of the "stars_data" columns columns which differ from the
  ;                     default values listed above. For example, column_names = {RA:'ICRS_RA'} can be used to specify that
  ;                     the RA values are listed in the column of stars_data named ICRS_RA.
  ;       RA - Right ascension (decimal degrees). A N-elements array can be specified to calculate the Bayesian probability of several
  ;            stars at once, but then all mandatory inputs must also be N-elements arrays.
  ;       DEC - Declination (decimal degrees).
  ;       PMRA - Proper motion in the right ascension direction (mas/yr, must include the cos(dec) factor).
  ;       PMDEC - Proper motion in the declination direction (mas/yr).
  ;       EPMRA - Measurement error on the proper motion in the right ascension direction (mas/yr, must not include the cos(dec) factor).
  ;       EPMDEC -  Measurement error on the proper motion in the declination direction (mas/yr).
  ;       RV - Radial velocity measurement to be included in the Bayesian probability (km/s).
  ;            If this keyword is set, ERV must also be set.
  ;            A N-elements array must be used if N stars are analyzed at once.
  ;       ERV - Measurement error on the radial velocity to be included in the Bayesian probability (km/s).
  ;             A N-elements array must be used if N stars are analyzed at once.
  ;       DIST - Distance measurement to be included in the Bayesian probability (pc).
  ;              By default, the BANYAN_SIGMA Bayesian priors are meant for this keyword to be used with trigonometric distances only.
  ;              Otherwise, the rate of true positives may be far from the nominal values described in Gagné et al. (ApJS, 2017, XX, XX).
  ;              If this keyword is set, EDIST must also be set.
  ;              A N-elements array must be used if N stars are analyzed at once.
  ;       EDIST - Measurement error on the distance to be included in the Bayesian probability (pc).
  ;               A N-elements array must be used if N stars are analyzed at once.
  ;       PLX - Parallax measurement to be included in the Bayesian probability (mas). The distance will be approximated with DIST = 1000/PLX.
  ;             If this keyword is set, EPLX must also be set.
  ;             A N-elements array must be used if N stars are analyzed at once.
  ;       EPLX - Measurement error on the parallax to be included in the Bayesian probability (mas). The distance error will be approximated with EDIST = 1000/PLX^2*EPLX.
  ;              A N-elements array must be used if N stars are analyzed at once.
  ;       PSIRA - Parallax motion factor PSIRA described in Gagné et al. (ApJS, 2017, XX, XX), in units of 1/yr.
  ;                If this keyword is set, the corresponding PSIDEC, EPSIRA and EPSIDEC keywords must also be set.
  ;                This measurement is only useful when proper motions are estimated from two single-epoch astrometric
  ;                measurements. It captures the dependence of parallax motion as a function of distance, and allows
  ;                BANYAN_SIGMA to shift the UVW center of the moving group models, which is equivalent to
  ;                correctly treating the input "proper motion" PMRA, PMDEC, EPMRA, EPMDEC as a true apparent motion.
  ;                This keyword should *not* be used if proper motions were derived from more than two epochs, or if
  ;                they were obtained from a full parallax solution.
  ;                A N-elements array must be used if N stars are analyzed at once.
  ;       PSIDEC - Parallax motion factor PSIDEC described in Gagné et al. (ApJS, 2017, XX, XX), in units of 1/yr.
  ;                 A N-elements array must be used if N stars are analyzed at once.
  ;       EPSIRA - Measurement error on the parallax motion factor PSIRA described in Gagné et al. (ApJS, 2017, XX, XX),
  ;                 in units of 1/yr. A N-elements array must be used if N stars are analyzed at once.
  ;       EPSIDEC - Measurement error on the parallax motion factor PSIDEC described in Gagné et al. (ApJS, 2017, XX, XX),
  ;                  in units of 1/yr. A N-elements array must be used if N stars are analyzed at once.
  ;       NTARGETS_MAX - (default 10^6). Maximum number of objects to run at once in BANYAN_SIGMA to avoid saturating the RAM.
  ;                      If more targets are supplied, BANYAN_SIGMA is run over a loop of several batches of NTARGETS_MAX objects. 
  ;       HYPOTHESES - The list of Bayesian hypotheses to be considered. They must all be present in the parameters fits file
  ;                    (See REQUIREMENTS #1 above).
  ;       LN_PRIORS - An IDL structure that contains the natural logarithm of Bayesian priors that should be *multiplied with the
  ;                   default priors* (use /UNIT_PRIORS if you want only LN_PRIORS to be considered). The structure must contain the name
  ;                   of each hypothesis as tags, and the associated scalar value of the natural logarithm of the Bayesian prior for each tag. 
  ;       CONSTRAINT_DIST_PER_HYP - An IDL structure (or array of IDL structures when several objects are analyzed) that contains a distance constraint (in pc).
  ;                   Each of the Bayesian hypotheses must be included as structure tags and the distance must be specified as its
  ;                   associated scalar value. CONSTRAINT_EDIST_PER_HYP must also be specified if CONSTRAINT_DIST_PER_HYP is specified.
  ;                   This keyword is useful for including spectro-photometric distance constraints that depend on the age of the young association or field.
  ;       CONSTRAINT_EDIST_PER_HYP - An IDL structure (or array of IDL structures when several objects are analyzed) that contains a measurement
  ;                   error on the distance constraint (in pc). Each of the Bayesian hypotheses must be included as structure tags and the
  ;                   distance error must be specified as its associated scalar value.  
  ;
  ; OPTIONAL INPUT KEYWORD:
  ;       /UNIT_PRIORS - If this keyword is set, all default priors are set to 1 (but they are still overrided by manual priors input with the keyword LN_PRIORS).
  ;       /LNP_ONLY - If this keyword is set, only Bayesian probabilities will be calculated and returned.
  ;       /NO_XYZ - If this keyword is set, the width of the spatial components of the multivariate Gaussian will be widened by a large
  ;                 factor, so that the XYZ components are effectively ignored. This keyword must be used with extreme caution as it will
  ;                 generate a significant number of false-positives and confusion between the young associations.
  ;       /USE_RV - Use any radial velocity values found in the stars_data input structure.
  ;       /USE_DIST - Use any distance values found in the stars_data input structure.
  ;       /USE_PLX - Use any parallax values found in the stars_data input structure.
  ;       /USE_PSI - Use any psira, psidec values found in the stars_data input structure.
  ;       /OVERRIDE_ERRORS - Do not exit program even when errors are encountered.
  ;
  ; OUTPUT:
  ;      This routine outputs a single IDL structure (or array of structures when many objects are analyzed at once), with the following tags:
  ;      NAME - The name of the object (as taken from the input structure).
  ;      ALL - A structure that contains the Bayesian probability (0 to 1) for each of the associations (as individual tags).
  ;      METRICS - A structure that contains the performance metrics associated with the global Bayesian probability of this target.
  ;                This sub-structure contains the following tags:
  ;        TPR - Rate of true positives expected in a sample of objects that have a Bayesian membership probability at least as large as that of the target.
  ;        FPR - Rate of false positives (from the field) expected in a sample of objects that have a Bayesian membership probability at least as large as that of the target.
  ;        PPV - Positive Predictive Value (sample contamination) expected in a sample of objects that have a Bayesian membership probability at least as large as that of the target.
  ;      [ASSOCIATION_1] - Sub-structure containing the relevant details for assiciation [ASSOCIATION_1].
  ;      [ASSOCIATION_2] - (...).
  ;      (...).
  ;      [ASSOCIATION_N] - (...).
  ;                        These sub-structures contain the following tags:
  ;        HYPOTHESIS - Name of the association
  ;        PROB - Bayesian probability (0 to 1)
  ;        D_OPT - Optimal distance (pc) that maximizes the Bayesian likelihood for this hypothesis.
  ;        RV_OPT - Optimal radial velocity (km/s) that maximizes the Bayesian likelihood for this hypothesis.
  ;        ED_OPT - Error on the optimal distance (pc), which approximates the 68% width of how the likelihood varies with distance.
  ;        ERV_OPT - Error on the optimal radial velocity (km/s), which approximates the 68% width of how the likelihood varies with
  ;                  radial velocity.
  ;        XYZUVW - 6-dimensional array containing the XYZ and UVW position of the star at the measured radial velocity and/or
  ;                 distance, or the optimal radial velocity and/or distance when the first are not available (units of pc and km/s).
  ;        EXYZUVW - Errors on XYZUVW (units of pc and km/s).
  ;        XYZ_SEP - Separation between the optimal or measured XYZ position of the star and the center of the multivariate Gaussian
  ;                  model of this Bayesian hypothesis (pc).
  ;        UVW_SEP - Separation between the optimal or measured UVW position of the star and the center of the multivariate Gaussian
  ;                  model of this Bayesian hypothesis (km/s).
  ;        XYZ_SEP - N-sigma separation between the optimal or measured XYZ position of the star and the multivariate Gaussian model
  ;                  of this Bayesian hypothesis (no units).
  ;        UVW_SEP - N-sigma separation between the optimal or measured UVW position of the star and the multivariate Gaussian model
  ;                  of this Bayesian hypothesis (no units).
  ;        MAHALANOBIS - Mahalanobis distance between the optimal or measured XYZUVW position of the star and the multivariate Gaussian
  ;                      model. A Mahalanobis distance is a generalization of a 6D N-sigma distance that accounts for covariances. 
  ;      BESTYA_STR - A sub-structure similar to those described above for the most probable young association (ignoring the field possibility).
  ;      YA_PROB - The Bayesian probability (0 to 1) that this object belongs to any young association (i.e., excluding the field).
  ;      LIST_PROB_YAS - A list of young associations with at least 5% Bayesian probability. Their relative probabilities (%) are specified
  ;                      between parentheses.
  ;      BEST_HYP - Most probable Bayesian hypothesis (including the field)
  ;      BEST_YA - Most probable single young association.
  ;      
  ; PROCEDURES USED:
  ;       BANYAN_SIGMA_SOLVE_MULTIVAR, NAN_STR, ALOG_SUM_2D, MRDFITS, REMOVE, UNIQ_UNSORTED, ADD_TAGS
  ;
  ; MODIFICATION HISTORY:
  ;       WRITTEN, Jonathan Gagne, October, 25 2017
  ;-
  
  forward_function banyan_sigma_solve_multivar, nan_str, alog_sum_2d, mrdfits, remove, uniq_unsorted, add_tags
  
  ;Maximal numer of targets to be put at once in the RAM for the analytical solving of the BANYAN-SIGMA integrals
  ; Using a number that is too high will cause IDL to crash
  if ~keyword_set(ntargets_max) then $
    ntargets_max = 1d6
  
  ;The total number of stars in the Besancon model within 300 pc to tranlate FPR to NFP
  total_besancon_objects = double(7152397L)
  
  if ~keyword_set(stars_data) and (ra eq !NULL or dec eq !NULL or pmra eq !NULL or pmdec eq !NULL or epmdec eq !NULL) then $
    message, ' Either an input structure (stars_data) or all of the ra,dec,pmra,pmdec,epmra and epmdec keywords must be specified !', continue=keyword_set(override_errors)
  
  if keyword_set(constraint_dist_per_hyp) and ~keyword_set(constraint_edist_per_hyp) then $
    message, ' If CONSTRAINT_DIST_PER_HYP is specified, CONSTRAINT_EDIST_PER_HYP must also be specified !', continue=keyword_set(override_errors)
  
  ;Default column names
  default_column_names = {RA:'RA',DEC:'DEC',PMRA:'PMRA',PMDEC:'PMDEC',$
    EPMRA:'EPMRA',EPMDEC:'EPMDEC'}
  if keyword_set(use_rv) then $
    add_tags, temporary(default_column_names), ['RV','ERV'], ['''RV''','''ERV'''], default_column_names
  if keyword_set(use_dist) then $
    add_tags, temporary(default_column_names), ['DIST','EDIST'], ['''DIST''','''EDIST'''], default_column_names
  if keyword_set(use_plx) then $
    add_tags, temporary(default_column_names), ['PLX','EPLX'], ['''PLX''','''EPLX'''], default_column_names
  if keyword_set(use_psi) then $
    add_tags, temporary(default_column_names), ['PSIRA','PSIDEC','EPSIRA','EPSIDEC'], ['''PSIRA''','''PSIDEC''','''EPSIRA''','''EPSIDEC'''], default_column_names
  
  ;Merge the default structure describing column names with the user-specified one
  if keyword_set(column_names) then begin
    user_tags = tag_names(column_names)
    ;default_tags = tag_names(default_column_names)
    for i=0L, n_elements(user_tags)-1L do begin
      gfound = where(tag_names(default_column_names) eq user_tags[i], ngfound)
      if ngfound eq 0L then $
        add_tags, temporary(default_column_names), user_tags[i], '''NaN''', default_column_names
      gtoadd = where(tag_names(default_column_names) eq user_tags[i], ngtoadd)
      default_column_names.(gtoadd[0L]) = column_names.(i)
    endfor
  endif
  column_names = default_column_names
  
  ;Check if a column named PLX, DIST, RV, PSIRA, etc. exist in stars_data but not in
  ;column_names. If this is the case, issue a warning so that the user understands that
  ;some data are not being considered.
  if keyword_set(stars_data) then begin
    if max(strpos(strupcase(tag_names(stars_data)),'RV') ne -1L) eq 1 and $
      max(strpos(strupcase(tag_names(column_names)),'RV') ne -1L) eq 0 and use_rv eq !NULL then $
      message, ' Warning: Radial velocities (RV) were not read from the input data, because the RV tag was not included in the COLUMN_NAMES keyword of BANYAN_SIGMA. You can also call BANYAN_SIGMA with the /USE_RV keyword to read them, or with USE_RV=0 to avoid this warning.', /continue
    if max(strpos(strupcase(tag_names(stars_data)),'DIST') ne -1L) eq 1 and $
      max(strpos(strupcase(tag_names(column_names)),'DIST') ne -1L) eq 0 and use_dist eq !NULL then $
      message, ' Warning: Distances (DIST) were not read from the input data, because the DIST tag was not included in the COLUMN_NAMES keyword of BANYAN_SIGMA. You can also call BANYAN_SIGMA with the /USE_DIST keyword to read them, or with USE_DIST=0 to avoid this warning.', /continue
    if max(strpos(strupcase(tag_names(stars_data)),'PLX') ne -1L) eq 1 and $
      max(strpos(strupcase(tag_names(column_names)),'PLX') ne -1L) eq 0 and use_plx eq !NULL then $
      message, ' Warning: Parallaxes (PLX) were not read from the input data, because the PLX tag was not included in the COLUMN_NAMES keyword of BANYAN_SIGMA. You can also call BANYAN_SIGMA with the /USE_PLX keyword to read them, or with USE_PLX=0 to avoid this warning.', /continue
    if (max(strpos(strupcase(tag_names(stars_data)),'PSIRA') ne -1L) eq 1 and $
      max(strpos(strupcase(tag_names(column_names)),'PSIRA') ne -1L) eq 0) or $
      (max(strpos(strupcase(tag_names(stars_data)),'PSIDEC') ne -1L) eq 1 and $
      max(strpos(strupcase(tag_names(column_names)),'PSIDEC') ne -1L) eq 0) and use_psi eq !NULL then $
      message, ' Warning: The PSI parameters (PSIRA,PSIDEC) were not read from the input data, because the PSIRA and PSIDEC tags were not included in the COLUMN_NAMES keyword of BANYAN_SIGMA. You can also call BANYAN_SIGMA with the /USE_PSI keyword to read them, or with USE_PSI=0 to avoid this warning.', /continue
  endif
  
  ;Read the input structure and avoid changing it
  s = {name:'NaN',RA:!values.d_nan,DEC:!values.d_nan,PMRA:!values.f_nan,PMDEC:!values.f_nan,$
    EPMRA:!values.f_nan,EPMDEC:!values.f_nan,RV:!values.f_nan,ERV:!values.f_nan,DIST:!values.f_nan,$
    EDIST:!values.f_nan,PSIRA:0d0,PSIDEC:0d0,EPSIRA:0d0,EPSIDEC:0d0}
  if keyword_set(stars_data) then begin
    nobj = n_elements(stars_data)
    s = replicate(s, nobj)
    
    ;Assign information from stars_data in the correct s entries
    ;stop
    column_tags = tag_names(column_names)
    for i=0L, n_tags(column_names)-1L do begin
      gsi = where(tag_names(s) eq column_tags[i], ngsi)
      if ngsi eq 0L then continue
      ;!!!!!!!!! JG NOV 6 2017 THE STRUPCASE WAS MISSING BELOW, MAKE SURE THIS IS FIXED IN PYTHON TOO !!!!!!!!!!!!
      gstarsi = where(tag_names(stars_data) eq strupcase(column_names.(i)), ngstarsi)
      if ngstarsi eq 0L then continue
      s.(gsi[0L]) = stars_data.(gstarsi[0L])
    endfor
  endif else begin
    nobj = n_elements(ra)
    s = replicate(s, nobj)
    s.RA = ra
    s.DEC = dec
    s.PMRA = pmra
    s.PMDEC = pmdec
    s.EPMRA = epmra
    s.EPMDEC = epmdec
  endelse
  
  ;Transform parallaxes to distances directly in the s structure
  if max(strpos(strupcase(tag_names(s)),'EPLX') ne -1L) eq 1 and max(strpos(strupcase(tag_names(s)),'PLX') ne -1L) eq 1 then begin
    if max(strpos(strupcase(tag_names(s)),'EDIST') ne -1L) eq 0 then $
      str_add_tags, temporary(s), ['EDIST'], ['!values.d_nan'], s, continue=keyword_set(override_errors)
    gfinite = where(finite(s.PLX) and finite(s.EPLX), ngfinite)
    if ngfinite ne 0L then $
      s[gfinite].EDIST = 1d3/s[gfinite].PLX^2*s[gfinite].EPLX
    s.EPLX = !values.d_nan
  endif
  if max(strpos(strupcase(tag_names(s)),'PLX') ne -1L) eq 1 then begin
    if max(strpos(strupcase(tag_names(s)),'DIST') ne -1L) eq 0 then $
      str_add_tags, temporary(s), ['DIST'], ['!values.d_nan'], s
    gfinite = where(finite(s.PLX), ngfinite)
    if ngfinite ne 0L then $
      s[gfinite].DIST = 1d3/s[gfinite].PLX
    s.PLX = !values.d_nan
  endif
  
  ;If measurements are specified as keywords, include them in the structure
  if ra ne !NULL then $
    s.ra = ra
  if dec ne !NULL then $
    s.dec = dec
  if pmra ne !NULL then $
    s.pmra = pmra
  if pmdec ne !NULL then $
    s.pmdec = pmdec
  if epmra ne !NULL then $
    s.epmra = epmra
  if epmdec ne !NULL then $
    s.epmdec = epmdec
  if dist ne !NULL then $
    s.dist = dist
  if edist ne !NULL then $
    s.edist = edist
  if plx ne !NULL then $
    s.dist = 1d3/plx
  if eplx ne !NULL and plx ne !NULL then $
    s.edist = 1d3/plx^2*eplx
  if rv ne !NULL then $
    s.rv = rv
  if erv ne !NULL then $
    s.erv = erv
  if psira ne !NULL then $
    s.psira = psira
  if psidec ne !NULL then $
    s.psidec = psidec
  if epsira ne !NULL then $
    s.epsira = epsira
  if epsidec ne !NULL then $
    s.epsidec = epsidec
  
  ;Check for unphysical data
  bad = where(s.ra lt 0 or s.ra ge 360., nbad)
  if nbad ne 0L then $
    message, ' Some RA values are unphysical !', continue=keyword_set(override_errors)
  
  bad = where(s.dec lt -90. or s.dec gt 90., nbad)
  if nbad ne 0L then $
    message, ' Some DEC values are unphysical !', continue=keyword_set(override_errors)
  
  bad = where(s.epmra le 0. or s.epmdec le 0., nbad)
  if nbad ne 0L then $
    message, ' Some EPMRA or EPMDEC values are unphysical !', continue=keyword_set(override_errors)
    
  bad = where(~finite(s.ra) or ~finite(s.dec) or ~finite(s.pmra) or ~finite(s.pmdec) or ~finite(s.epmra) or ~finite(s.epmdec), nbad)
  if nbad ne 0L then $
    message, ' The observables ra,dec,pmra,pmdec,epmra and epmdec must be specified (and finite) for each object !', continue=keyword_set(override_errors)
  
  if max(strpos(strupcase(tag_names(s)),'ERV') ne -1L) eq 1 then begin
    bad = where(s.erv le 0., nbad)
    if nbad ne 0L then $
      message, ' Some ERV values are unphysical !', continue=keyword_set(override_errors)
  endif
  
  if max(strpos(strupcase(tag_names(s)),'DIST') ne -1L) eq 1 and $
    max(strpos(strupcase(tag_names(s)),'EDIST') ne -1L) eq 1 then begin
    bad = where(s.dist lt 0. or s.edist le 0., nbad)
    if nbad ne 0L then $
      message, ' Some DIST or EDIST values are unphysical !', continue=keyword_set(override_errors)
  endif
  
  bad = where(((s.psira ne 0. or s.psidec ne 0.) and (s.epsira eq 0. or s.epsidec eq 0.)) or s.epsira lt 0. or s.epsidec lt 0., nbad)
  if nbad ne 0L then $
    message, ' Some EPSIRA or EPSIDEC values are unphysical !', continue=keyword_set(override_errors)
  
  bad = where(finite(s.dist) and ~finite(s.edist), nbad)
  if nbad ne 0L then $
    message, ' Some DIST values are specified without EDIST !', continue=keyword_set(override_errors)
  
  bad = where(finite(s.rv) and ~finite(s.erv), nbad)
  if nbad ne 0L then $
    message, ' Some RV values are specified without ERV !', continue=keyword_set(override_errors)
  
  ;Locate the parameters file
  parameters_file = file_dirname(routine_filepath())+path_sep()+'data'+path_sep()+'banyan_sigma_parameters.fits'
  if ~file_test(parameters_file) then $
    message, ' The multivariate Gaussian parameters file could not be found ! Please make sure that you did not move "'+path_sep()+'data'+path_sep()+'banyan_sigma_parameters.fits" from the same path as the IDL routine banyan_sigma.pro !', continue=keyword_set(override_errors)
  
  ;Read the multivariate Gaussian parameters describing the spatial-kinematic structure of all hypotheses
  parameters_str = mrdfits(parameters_file,1,/silent)
  parameters_str.name = strtrim(parameters_str.name,2)
  npar = n_elements(parameters_str)
  
  ;Replace parameter names that start with a number to avoid problems with using them as structure tags
  bad = where(is_number(strmid(parameters_str.NAME,0,1)), nbad)
  if nbad ne 0L then parameters_str[bad].NAME = '_'+parameters_str[bad].NAME
  
  ;Build a list of the default hypotheses to be considered
  if ~keyword_set(hypotheses) then $
    hypotheses = strtrim(parameters_str[uniq_unsorted(parameters_str.name)].name,2)
  nhyp = n_elements(hypotheses)

  ;Replace hypothesis names that start with a number to avoid problems with using them as structure tags
  bad = where(is_number(strmid(hypotheses,0,1)), nbad)
  if nbad ne 0L then hypotheses[bad] = '_'+hypotheses[bad]
  
  ;If CONSTRAINT_DIST_PER_HYP is set, check that all hypotheses are included
  if keyword_set(constraint_dist_per_hyp) then begin
    dist_per_hyp_tags = tag_names(constraint_dist_per_hyp)
    edist_per_hyp_tags = tag_names(constraint_edist_per_hyp)
    if ~array_equal(dist_per_hyp_tags[sort(dist_per_hyp_tags)],hypotheses[sort(hypotheses)]) or ~array_equal(edist_per_hyp_tags[sort(dist_per_hyp_tags)],hypotheses[sort(hypotheses)]) then begin
      print, ' The tag names of CONSTRAINT_DIST_PER_HYP are not equal to the set of hypotheses :'
      print, transpose('  '+hypotheses)
      message, ' Please correct the tag names appropriately', continue=keyword_set(override_errors)
    endif
    
    ;Build CONSTRAINT_DIST_PER_HYP into an array
    dist_per_hyp_arr = make_array(nobj,nhyp,value=!values.d_nan,/double)
    edist_per_hyp_arr = make_array(nobj,nhyp,value=!values.d_nan,/double)
    for i=0L, nhyp-1L do begin
      gtag = where(dist_per_hyp_tags eq hypotheses[i] or dist_per_hyp_tags eq '_'+hypotheses[i], ngtag)
      dist_per_hyp_arr[*,i] = double(constraint_dist_per_hyp.(gprior[0L]))
    endfor
    for i=0L, nhyp-1L do begin
      gtag = where(edist_per_hyp_tags eq hypotheses[i] or edist_per_hyp_tags eq '_'+hypotheses[i], ngtag)
      edist_per_hyp_arr[*,i] = double(constraint_edist_per_hyp.(gprior[0L]))
    endfor
    
    ;Check that all distance constraints are physical
    bad = where(dist_per_hyp_arr lt 0. or edist_per_hyp_arr le 0., nbad)
    if nbad ne 0L then $
      message, ' Some of the specified CONSTRAINT_DIST_PER_HYP or CONSTRAINT_EDIST_PER_HYP values are unphysical !', continue=keyword_set(override_errors)
    
    ;Check that the CONSTRAINT_EDIST_PER_HYP values are finite everywhere that CONSTRAINT_DIST_PER_HYP are finite
    bad = where(finite(dist_per_hyp_arr) and ~finite(edist_per_hyp_arr), nbad)
    if nbad ne 0L then $
      message, ' Some of the specified CONSTRAINT_EDIST_PER_HYP are not finite where CONSTRAINT_DIST_PER_HYP are finite !', continue=keyword_set(override_errors)
    
    ;Check that either all or none of the distance constraints are finite for a given object
    bad = where((finite(total(dist_per_hyp_arr,2,/nan)) and ~finite(total(dist_per_hyp_arr,2))) or (finite(total(edist_per_hyp_arr,2,/nan)) and ~finite(total(edist_per_hyp_arr,2))), nbad)
    if nbad ne 0L then $
      message, ' The CONSTRAINT_DIST_PER_HYP and CONSTRAINT_EDIST_PER_HYP values must be all finite or all non-finite for a given star !', continue=keyword_set(override_errors)
  endif
  
  ;Override priors to all ones if UNIT_PRIORS is set
  if keyword_set(unit_priors) then $
    parameters_str.ln_prior = 0d0
  
  ;Determine whether a trigonometric distance or a per-hypothesis distance constraint was set 
  if keyword_set(constraint_dist_per_hyp) then $
    distance_is_set = finite(s.dist) or finite(total(dist_per_hyp_arr,2,/nan)) else $
    distance_is_set = finite(s.dist)
  
  ;Adjust the priors
  g_pm = where(~finite(s.rv) and ~distance_is_set, ng_pm)
  g_pm_rv = where(finite(s.rv) and ~distance_is_set, ng_pm_rv)
  g_pm_dist = where(~finite(s.rv) and distance_is_set, ng_pm_dist)
  g_pm_rv_dist = where(finite(s.rv) and distance_is_set, ng_pm_rv_dist)
  ln_priors_nd = make_array(nobj,nhyp,value=0d0,/double)
  ln_priors_nd_manual = make_array(nobj,nhyp,value=0d0,/double)
  for i=0L, nhyp-1L do begin
    ;Skip the field hypotheses as they do not have a Bayesian prior
    if strpos(strupcase(hypotheses[i]),'FIELD') ne -1L then continue
    ;Read the parameters structure to identify the 4 priors associated with a given young association
    ind = where(strtrim(strlowcase(parameters_str.name),2) eq strtrim(strlowcase(hypotheses[i]),2) or '_'+strtrim(strlowcase(parameters_str.name),2) eq strtrim(strlowcase(hypotheses[i]),2), nfi)
    if nfi ne 1 then $
      message, ' Hypothesis '+hypotheses[i]+' could not be found in the parameters structure !'
    ln_priors_i = parameters_str[ind[0L]].LN_PRIOR
    ;In the cases where only one prior is designated, assign it to all stars
    if n_elements(ln_priors_i) eq 1L then $
      ln_priors_nd[*,i] = ln_priors_i[0L]
    ;Otherwise assign them properly as a function of available observables
    if n_elements(ln_priors_i) eq 4L then begin
      if ng_pm ne 0L then $
        ln_priors_nd[g_pm,i] = ln_priors_i[0L]
      if ng_pm_rv ne 0L then $
        ln_priors_nd[g_pm_rv,i] = ln_priors_i[1L]
      if ng_pm_dist ne 0L then $
        ln_priors_nd[g_pm_dist,i] = ln_priors_i[2L]
      if ng_pm_rv_dist ne 0L then $
        ln_priors_nd[g_pm_rv_dist,i] = ln_priors_i[3L]
    endif
  endfor
  
  ;Add the manual priors if they are specified as an input structure
  if keyword_set(ln_priors) then begin
    for i=0L, nhyp-1L do begin
      ;The field hypotheses *can* have manual priors
      gprior = where(tag_names(ln_priors) eq hypotheses[i] or tag_names(ln_priors) eq '_'+hypotheses[i], ngprior)
      if ngprior eq 0 then begin
        message, ' The prior for hypothesis '+hypotheses[i]+' was left to its default value as it was not specified manually', /continue
        continue
      endif
      ln_priors_nd_manual[*,i] = double(ln_priors.(gprior[0L]))
    endfor
    
    ;Normalize manual priors with the field hypothesis (because they get applied only on young associations)
    gnorm = where(strpos(strupcase(hypotheses),'FIELD') ne -1L, ngnorm)
    if ngnorm eq 1 then $
      norm_priors_1d = ln_priors_nd_manual[*,gnorm[0L]] else $
      norm_priors_1d = alog_sum_2d(ln_priors_nd_manual[*,gnorm],dim=2)
    ln_priors_nd_manual -= (norm_priors_1d#make_array(nhyp,value=1d0,/double))
    
    ;Apply the manual priors on top of the default priors
    ln_priors_nd += ln_priors_nd_manual
  endif
  
  ;If both trigonometric distances and per-hypothesis distance constraints are set, transform the per-hypothesis distance constraints into priors
  nboth_distances_set = 0L
  if keyword_set(constraint_dist_per_hyp) then $
    both_distances_set = where(finite(s.dist) and finite(total(dist_per_hyp_arr,2,/nan)), nboth_distances_set)
  if nboth_distances_set ne 0L then begin
    ones_hyp = make_array(nhyp,value=1d0,/double)
    ln_prob_dist_differences = -((s[both_distances_set].dist#ones_hyp)-dist_per_hyp_arr[both_distances_set,*])^2/(2d0*(((s[both_distance_set].edist^2d0)#ones_hyp)+edist_per_hyp_arr[both_distances_set,*]^2))
    
    ;Treat these values as priors so normalize them with the field hypotheses (because they get applied only on young associations)
    gnorm = where(strpos(strupcase(hypotheses),'FIELD') ne -1L, ngnorm)
    if ngnorm eq 1 then $
      norm_priors_1d = ln_prob_dist_differences[*,gnorm[0L]] else $
      norm_priors_1d = alog_sum_2d(ln_prob_dist_differences[*,gnorm],dim=2)
    ln_prob_dist_differences -= (norm_priors_1d#ones_hyp)
    
    ;Apply these values on the priors
    ln_priors_nd[both_distances_set,*] += ln_prob_dist_differences
    
    ;Remove the per-hypothesis distance constraints on these particular objects and just keep the trigonometric distances
     dist_per_hyp_arr[both_distances_set,*] = !values.d_nan
     edist_per_hyp_arr[both_distances_set,*] = !values.d_nan
  endif
  
  ;Initiate an array that will contain the ln probabilities if those are the only required outputs  
  if keyword_set(lnp_only) then $
    all_lnprobs = dblarr(nobj,nhyp)*!values.d_nan
  
  ;Loop on hypotheses to run BANYAN Sigma on
  for i=0L, nhyp-1L do begin
    ;Find the right spatial-kinematic model and store it in "association_structure"
    ind = where(strtrim(strlowcase(parameters_str.name),2) eq strtrim(strlowcase(hypotheses[i]),2) or '_'+strtrim(strlowcase(parameters_str.name),2) eq strtrim(strlowcase(hypotheses[i]),2), ngauss)
    ;ngauss indicates the number of multivariate Gaussian parameters
    
    ;If CONSTRAINT_DIST_PER_HYP is set, determine which distance constraint must be used now
    dist_for_this_hypothesis = s.dist
    edist_for_this_hypothesis = s.edist
    if keyword_set(constraint_dist_per_hyp) then begin
      gdist_per_hyp = where(finite(constraint_dist_per_hyp[*,i]), ngdist_per_hyp)
      if ngdist_per_hyp ne 0L then begin
        dist_for_this_hypothesis = constraint_dist_per_hyp[gdist_per_hyp,i]
        edist_for_this_hypothesis = constraint_edist_per_hyp[gdist_per_hyp,i]
      endif
    endif
    
    ;Loop on individual multivariate gaussians, if needed
    output_str_multimodel = !NULL
    if keyword_set(lnp_only) then $
      all_lnprobs_hypi = dblarr(nobj,ngauss)
    for gaussi=0L, ngauss-1L do begin
      
      association_structure = parameters_str[ind[gaussi]]
      association_structure.name = strtrim(association_structure.name,2)
      
      ;Determine how many batches will be needed to avoid saturating the RAM
      nbatches = ceil(double(nobj)/double(ntargets_max))
      for ci=0L, nbatches-1L do begin

        ;Determine the indices of targets to be selected
        ind_from = round(ci * ntargets_max)
        ind_to = ind_from + round(ntargets_max-1L)
        ind_to <= round(nobj-1L)

        ;Create a sub-structure of input data
        s_ci = s[ind_from:ind_to]
        dist_for_this_hypothesis_ci = dist_for_this_hypothesis[ind_from:ind_to]
        edist_for_this_hypothesis_ci = edist_for_this_hypothesis[ind_from:ind_to]
        nobj_ci = n_elements(s_ci)
        
        
        ;Solve the BANYAN-SIGMA integrals for one hypothesis and N targets
        output_str_ci = banyan_sigma_solve_multivar(s_ci.ra, s_ci.dec, s_ci.pmra, s_ci.pmdec, s_ci.epmra, s_ci.epmdec, RV_MES=s_ci.rv, ERV_MES=s_ci.erv, $
          DIST_MES=dist_for_this_hypothesis_ci, EDIST_MES=edist_for_this_hypothesis_ci, ASSOCIATION_STRUCTURE=association_structure, PSIRA=s_ci.psira, $
          PSIDEC=s_ci.psidec, EPSIRA=s_ci.epsira, EPSIDEC=s_ci.epsidec, LNP_ONLY=lnp_only, NO_XYZ=no_xyz)
        
        if keyword_set(lnp_only) then begin
          all_lnprobs_hypi[ind_from:ind_to,gaussi] = output_str_ci.LN_P
          continue
        endif
        
        ;Create a general output structure if it does not exist
        if ~keyword_set(output_str) then begin
          output_str = output_str_ci[0L]
          nan_str, output_str
          output_str = replicate(output_str, nobj)
        endif

        ;Fill up the output structure
        output_str[ind_from:ind_to] = output_str_ci

      endfor;End of loop on chunks of data that avoids RAM saturation
      output_str_sci = !NULL
      s_ci = !NULL

      ;Skip the structure management if only a minimal output is required
      if keyword_set(lnp_only) then $
        continue

      ;Reformat the output structure if this hypothesis is a multivariate Gaussian mixture
      if ngauss ne 1 then begin
        if ~keyword_set(output_str_multimodel) then begin
          output_str_multimodel = output_str[0L]
          nan_str, output_str_multimodel
          output_str_multimodel = replicate(output_str_multimodel, nobj, ngauss)
        endif
        output_str_multimodel[*,gaussi] = output_str
        output_str = !NULL
      endif
    
    endfor;End of loop on individual multivariate gaussian components in the model
    
    ;If only log probs are required, compile them in the main array
    if keyword_set(lnp_only) then begin
      if ngauss eq 1L then $
        all_lnprobs[*,i] = all_lnprobs_hypi else $
      begin
        weights = parameters_str[ind].coefficient/total(parameters_str[ind].coefficient,/nan)
        alogweights2D = make_array(nobj,value=1d0,/double)#alog(weights)
        all_lnprobs[*,i] = alog_sum_2d(alogweights2D+all_lnprobs_hypi,dim=2L)
      endelse
      continue
    endif
    
    ;The output structure needs to be reformatted if there are more than one multivariate gaussians
    if ngauss ne 1 then begin
      output_str = output_str_multimodel[0L]
      nan_str, output_str
      output_str = replicate(output_str, nobj)
      weights = parameters_str[ind].coefficient/total(parameters_str[ind].coefficient,/nan)
      weights2D = make_array(nobj,value=1d0,/double)#weights
      output_str.LN_P = alog_sum_2d(alog(weights2D)+output_str_multimodel.LN_P,dim=2L)
      output_str.D_OPT = total(weights2D*output_str_multimodel.D_OPT,2L,/nan)
      output_str.ED_OPT = total(weights2D*output_str_multimodel.ED_OPT,2L,/nan)
      output_str.RV_OPT = total(weights2D*output_str_multimodel.RV_OPT,2L,/nan)
      output_str.ERV_OPT = total(weights2D*output_str_multimodel.ERV_OPT,2L,/nan)
      output_str.XYZ_SEP = total(weights2D*output_str_multimodel.XYZ_SEP,2L,/nan)
      output_str.UVW_SEP = total(weights2D*output_str_multimodel.UVW_SEP,2L,/nan)
      output_str.XYZ_SIG = total(weights2D*output_str_multimodel.XYZ_SIG,2L,/nan)
      output_str.UVW_SIG = total(weights2D*output_str_multimodel.UVW_SIG,2L,/nan)
      output_str.MAHALANOBIS = total(weights2D*output_str_multimodel.MAHALANOBIS,2L,/nan)
    endif
    output_str_multimodel = !NULL

    ;Create an output structure for all hypotheses and targets
    if ~keyword_set(output_str_all) then begin
      output_str_all = output_str[0L]
      nan_str, output_str_all
      output_str_all = replicate(output_str_all, nobj, nhyp)
    endif

    ;Store all
    output_str_all[*,i] = output_str
    output_str = !NULL
    
  endfor
  
  ;Useful shortcuts to calculate normalized probabilities
  ovec = make_array(nobj,value=1d0,/double)
  hvec = make_array(nhyp,value=1d0,/double)
  ovec = !NULL
  
  ;Store probabilities in a single array if needed
  if ~keyword_set(lnp_only) then $
    all_lnprobs = output_str_all.LN_P
  
  ;Normalize probabilities directly in log space
  ln_norm_output = all_lnprobs-(alog_sum_2d(all_lnprobs,dim=2L)#hvec)
  
  ;Compute [0,1] probabilities
  norm_output = exp(ln_norm_output)
  
  ;Identify hypotheses that correspond to moving groups or associations
  yind = where(strpos(strupcase(hypotheses),'FIELD') eq -1L, nyind, complement=ffind, ncomplement=nffind)
  
  ;Create an array of normalized YMG probabilities (no field)
  ln_norm_output_only_ymg = all_lnprobs[*,yind]-(alog_sum_2d(all_lnprobs[*,yind],dim=2L)#(hvec[yind]))

  ;Calculate the weighted YMG prior
  ln_prior_moving_groups = alog_sum_2d(ln_priors_nd[*,yind]+ln_norm_output_only_ymg,dim=2L)
  
  ;Weight the priors w/r/t the Bayesian probabilities and project these priors onto the field. This is a way to
  ; avoid having the priors change the relative moving group probabilities, as their goal is strictly to
  ; maximize young association vs FIELD classification performance
  ;Normalize probabilities directly in log space, projecting the inverse young association prior on the field probability
  ln_P_with_prior = all_lnprobs
  ln_P_with_prior[*,ffind] -= (ln_prior_moving_groups#make_array(nffind,value=1d0,/double))
  ;Renormalize
  ln_norm_output_prior = ln_P_with_prior-(alog_sum_2d(ln_P_with_prior,dim=2L)#hvec)
  
  ;Return log probabilities if this is the only required output
  if keyword_set(lnp_only) then $
    return, ln_norm_output_prior
  
  ;Compute [0,1] probabilities
  norm_output_prior = exp(ln_norm_output_prior)
  
  ;Read the association performance metrics
  ;Locate the parameters file
  metrics_file = file_dirname(routine_filepath())+path_sep()+'data'+path_sep()+'banyan_sigma_metrics.fits'
  metrics_computed = 0
  if ~file_test(metrics_file) then $
    message, ' The performance metrics file could not be found ! Performance metrics will not be calculated. Please make sure that you did not move "'+os.sep+'data'+os.sep+'banyan_sigma_metrics.fits" from the same path as the IDL file banyan_sigma.pro !', /continue
  ;Avoid computing biased metrics if the unit_priors keyword was set
  if file_test(metrics_file) and ~keyword_set(unit_priors) then begin 
    metrics_str = mrdfits(metrics_file,1,/silent)
    metrics_str.name = strtrim(metrics_str.name,2)
    
    ;Loop on young associations to determine their individual metrics
    tpr = dblarr(nobj,nyind)+!values.d_nan
    fpr = dblarr(nobj,nyind)+!values.d_nan
    ppv = dblarr(nobj,nyind)+!values.d_nan
    for yindi=0L, nyind-1L do begin
      gyindi = where(strtrim(metrics_str.name,2) eq hypotheses[yind[yindi]] or '_'+strtrim(metrics_str.name,2) eq hypotheses[yind[yindi]], ngyindi)
      if ngyindi gt 1 then $
        message, ' Bayesian hypothesis "'+hypotheses[yind[yindi]]+'" was found more than once in the performance metrics file ! Performance metrics should never be calculated for the FIELD hypothesis.'
      if ngyindi eq 0 then $
        message, ' Bayesian hypothesis "'+hypotheses[yind[yindi]]+'" was not found in the performance metrics file ! Its performance metrics will not be calculated. This will cancel out performance metrics calculations for any object with a >= 1% probability in that group.'
      
      probs_yindi = exp(ln_norm_output_prior[*,yindi] - alog_sum_2d([[ln_norm_output_prior[*,yindi]],[ln_norm_output_prior[*,ffind[0L]]]],dim=2))
      
      if ng_pm ne 0L then begin
        tpr[g_pm,yindi] = interpol(metrics_str[gyindi[0L]].tpr[*,0],metrics_str[gyindi[0L]].probs,probs_yindi[g_pm])
        fpr[g_pm,yindi] = interpol(metrics_str[gyindi[0L]].fpr[*,0],metrics_str[gyindi[0L]].probs,probs_yindi[g_pm])
        ppv[g_pm,yindi] = interpol(metrics_str[gyindi[0L]].ppv[*,0],metrics_str[gyindi[0L]].probs,probs_yindi[g_pm])
      endif
      if ng_pm_rv ne 0L then begin
        tpr[g_pm_rv,yindi] = interpol(metrics_str[gyindi[0L]].tpr[*,1],metrics_str[gyindi[0L]].probs,probs_yindi[g_pm_rv])
        fpr[g_pm_rv,yindi] = interpol(metrics_str[gyindi[0L]].fpr[*,1],metrics_str[gyindi[0L]].probs,probs_yindi[g_pm_rv])
        ppv[g_pm_rv,yindi] = interpol(metrics_str[gyindi[0L]].ppv[*,1],metrics_str[gyindi[0L]].probs,probs_yindi[g_pm_rv])
      endif
      if ng_pm_dist ne 0L then begin
        tpr[g_pm_dist,yindi] = interpol(metrics_str[gyindi[0L]].tpr[*,2],metrics_str[gyindi[0L]].probs,probs_yindi[g_pm_dist])
        fpr[g_pm_dist,yindi] = interpol(metrics_str[gyindi[0L]].fpr[*,2],metrics_str[gyindi[0L]].probs,probs_yindi[g_pm_dist])
        ppv[g_pm_dist,yindi] = interpol(metrics_str[gyindi[0L]].ppv[*,2],metrics_str[gyindi[0L]].probs,probs_yindi[g_pm_dist])
      endif
      if ng_pm_rv_dist ne 0L then begin
        tpr[g_pm_rv_dist,yindi] = interpol(metrics_str[gyindi[0L]].tpr[*,3],metrics_str[gyindi[0L]].probs,probs_yindi[g_pm_rv_dist])
        fpr[g_pm_rv_dist,yindi] = interpol(metrics_str[gyindi[0L]].fpr[*,3],metrics_str[gyindi[0L]].probs,probs_yindi[g_pm_rv_dist])
        ppv[g_pm_rv_dist,yindi] = interpol(metrics_str[gyindi[0L]].ppv[*,3],metrics_str[gyindi[0L]].probs,probs_yindi[g_pm_rv_dist])
      endif
    endfor
    
    ;Build the combination weights
    ln_weights = ln_norm_output_only_ymg
    ;Any group with less than 1% probability is ignored to avoid propagating potential NaNs
    bad = where(ln_weights lt alog(1d-2), nbad)
    ln_weights[bad] = !values.d_nan
    ;Re-normalize weights
    ln_weights -= (alog_sum_2d(ln_weights,dim=2L)#make_array(nyind,value=1d0,/double))
    
    ;Calculate the weighted metrics
    tpr_weighted = exp(alog_sum_2d(alog(tpr)+ln_weights,dim=2L))
    fpr_weighted = exp(alog_sum_2d(alog(fpr)+ln_weights,dim=2L))
    ppv_weighted = exp(alog_sum_2d(alog(ppv)+ln_weights,dim=2L))
    metrics_computed = 1
    
  endif
  
  ;Calculate the probability that each star is not part of the field
  hvec2 = make_array(nyind,value=1d0,/double)
  P_notfield = reform(exp(output_str_all[*,yind].LN_P),nobj,nyind)
  norm_output2 = P_notfield/(total(P_notfield,2,/nan)#hvec2)
  
  ;Determine the most probable hypothesis
  void = max(norm_output_prior,dim=2,/nan,wmax)
  most_probable_index = reform((array_indices(norm_output_prior,wmax))[1L,*])
  ;ghyp = where(strpos(strupcase(hypotheses),'FIELD') eq -1L, nghyp)
  if nyind gt 1 then $
    nya_prob_prior = total(norm_output_prior[*,yind],2,/nan)
  if nyind eq 1 then $
    nya_prob_prior = norm_output_prior[*,yind[0L]]
  
  ;Create a list of non-field hypotheses
  nn = !values.d_NAN
  hypotheses2 = strupcase(hypotheses[yind])
  ;Create an output structure that will receive all data
  void = execute('sout = replicate({NAME:''NaN'', ALL:{'+strjoin(strupcase(hypotheses)+':nn',', ')+'}, '+$
    'METRICS:{TPR:nn,FPR:nn,PPV:nn,NFP:nn}, '+$
    strjoin(strupcase([hypotheses,'BESTYA_STR'])+':{HYPOTHESIS:''NaN'',PROB:nn, D_OPT:nn, ED_OPT:nn, RV_OPT:nn, ERV_OPT:nn, '+$
    'XYZUVW:replicate(nn,6), EXYZUVW:replicate(nn,6), XYZ_SEP:nn, UVW_SEP:nn, XYZ_SIG:nn, UVW_SIG:nn, MAHALANOBIS:nn}',', ')+$
    ', YA_PROB:nn, LIST_PROB_YAS:''NaN'', BEST_HYP:''NaN'', BEST_YA:''NaN''}, nobj)')
  if void eq 0 then message, ' An execution statement has failed !'

  ;Fill up basic output info in the structure
  ;Name of target
  sout.NAME = s.NAME
  ;Probability that it belongs to any young association
  sout.YA_PROB = nya_prob_prior
  ;Most probable hypothesis
  sout.BEST_HYP = hypotheses[most_probable_index]
  
  ;Store performance metrics
  if metrics_computed eq 1 then begin
    sout.metrics.TPR = tpr_weighted
    sout.metrics.FPR = fpr_weighted
    sout.metrics.PPV = ppv_weighted
    sout.metrics.NFP = fpr_weighted*total_besancon_objects
  endif
  
  ;Fill up all associations sub-structures
  bbad = []
  for i=0L, nhyp-1L do begin
    h = hypotheses[i]
    ;Probability in list of all YMG probs
    void *= execute('sout.ALL.'+h+' = norm_output_prior[*,i]')
    ;Name of association
    void *= execute('sout.'+h+'.HYPOTHESIS = '''+h+'''')
    ;Probability of association
    void *= execute('sout.'+h+'.PROB = norm_output_prior[*,i]')
    ;Probability of association (priors=1)
    ;Statistical distance for this association
    void *= execute('sout.'+h+'.D_OPT = output_str_all[*,i].D_OPT')
    ;Error bars on statistical distance for this association
    void *= execute('sout.'+h+'.ED_OPT = output_str_all[*,i].ED_OPT')
    ;Statistical RV for this association
    void *= execute('sout.'+h+'.RV_OPT = output_str_all[*,i].RV_OPT')
    ;Error bars on statistical RV for this association
    void *= execute('sout.'+h+'.ERV_OPT = output_str_all[*,i].ERV_OPT')
    ;Statistical XYZ UVW position for this association
    void *= execute('sout.'+h+'.XYZUVW = output_str_all[*,i].XYZUVW')
    ;Error bars on statistical XYZ UVW position for this association
    void *= execute('sout.'+h+'.EXYZUVW = output_str_all[*,i].EXYZUVW')
    ;Separation from statistical XYZ position to the center of this association
    void *= execute('sout.'+h+'.XYZ_SEP = output_str_all[*,i].XYZ_SEP')
    ;Separation from statistical UVW position to the center of this association
    void *= execute('sout.'+h+'.UVW_SEP = output_str_all[*,i].UVW_SEP')
    ;Separation from statistical XYZ position to the center of this association divided by a quadrature sum of the spatial size of the association and the error bars on the statistical XYZ
    void *= execute('sout.'+h+'.XYZ_SIG = output_str_all[*,i].XYZ_SIG')
    ;Separation from statistical UVW position to the center of this association divided by a quadrature sum of the kinematic size of the association and the error bars on the statistical UVW
    void *= execute('sout.'+h+'.UVW_SIG = output_str_all[*,i].UVW_SIG')
    ;Mahalanobis (n-sigma for multivariate gaussians) distance between the statistical XYZUVW position and the spatial-kinematic model center of the association
    void *= execute('sout.'+h+'.MAHALANOBIS = output_str_all[*,i].MAHALANOBIS')
    if void eq 0 then message, ' An execution statement has failed !'
    nbad = 0L
    ;Flag objects with NaN probabilities
    void *= execute('bad = where(~finite(sout.'+h+'.PROB), nbad)')
    if void eq 0 then message, ' An execution statement has failed !'
    if nbad ne 0L then bbad = [bbad, bad]
  endfor
  norm_output = !NULL
  norm_output_prior = !NULL
  output_str_all = !NULL

  ;Loop on objects that is needed for tags that depend on the individual winning hypotheses
  for j=0L, nobj-1L do begin

    ;Create a string that lists all possible associations with more than 5% probability
    gg = where(reform(norm_output2[j,*]) gt .05, ngg)
    ;Sort the gg indices in order of decreasing probabilities
    gg = gg[reverse(sort(norm_output2[j,gg]))]
    if ngg eq 0L then continue
    if total(reform(norm_output2[j,*]),/nan) eq 0 then continue
    if ngg eq 1L then $
      sout[j].LIST_PROB_YAS = strupcase(hypotheses2[gg]) else $
      sout[j].LIST_PROB_YAS = strjoin(strupcase(hypotheses2[gg])+'('+strtrim(round(norm_output2[j,gg]*100.),2)+')',';')
    vmax = max(reform(norm_output2[j,*]),wmax,/nan)

    ;Store the most probable association structure in the BESTYA_STR tag
    void = execute('winh = sout[j].'+hypotheses2[wmax])
    sout[j].BEST_YA = hypotheses2[wmax]
    if void ne 0 then begin
      sout[j].BESTYA_STR = winh
    endif else begin
      message, ' An execution statement has failed in the determination of the WINH_STRUCTURE !', /continue
    endelse
  endfor
  norm_output2 = !NULL

  ;Fix the results of objects that have NaN probabilities
  if n_elements(bbad) ne 0L then begin
    sout[bbad].BEST_HYP = 'NaN'
    sout[bbad].YA_PROB = !values.f_nan
    for i=0L, nhyp-1 do begin
      void = execute('sout[bbad].'+hypotheses[i]+'.PROB = !values.f_nan')
      if void eq 0 then message, ' An execution statement has failed !'
    endfor
  endif
  
  ;Return the structure that contains all information on membership probabilities
  return, sout

End
