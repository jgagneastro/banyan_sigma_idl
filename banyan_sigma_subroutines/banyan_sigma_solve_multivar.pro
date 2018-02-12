Function banyan_sigma_solve_multivar, ra, dec, pmra, pmdec, epmra, epmdec, PRECISION_MATRIX=precision_matrix_in, $
  PRECISION_DETERM=precision_determ, CENTER_VECTOR=center_vector, RV_MES=rv_mes, ERV_MES=erv_mes, DIST_MES=dist_mes, $
  EDIST_MES=edist_mes, ASSOCIATION_STRUCTURE=association_structure, PSIRA=psira, PSIDEC=psidec, EPSIRA=epsira, $
  EPSIDEC=epsidec, LNP_ONLY=lnp_only, NO_XYZ=no_xyz

  ;+
  ; NAME:
  ;       BANYAN_SIGMA_SOLVE_MULTIVAR
  ;
  ; PURPOSE:
  ;       Solve the radial velocity and distance marginalization integrals (if needed) and compute log(probability) with
  ;       Bayes theorem for an array of stars and a single multivariate Gaussian XYZUVW model. This is a subroutine of banyan_sigma.pro
  ;       
  ;       Please see Gagne et al. 2017 (ApJS, XX, XX; Arxiv YY, YY) for more detail.
  ;       
  ;       An online version of this tool is available for 1-object queries at http://www.astro.umontreal.ca/~gagne/banyansigma.php.
  ;
  ; CALLING SEQUENCE:
  ;       OUTPUT_STRUCTURE = BANYAN_SIGMA_SOLVE_MULTIVAR(ra, dec, pmra, pmdec, epmra, epmdec, 
  ;         [ PRECISION_MATRIX=precision_matrix, PRECISION_DETERM=precision_determ, CENTER_VECTOR=center_vector,
  ;           RV_MES=rv_mes, ERV_MES=erv_mes, DIST_MES=dist_mes, EDIST_MES=edist_mes, 
  ;           ASSOCIATION_STRUCTURE=association_structure, PSIRA=psira, PSIDEC=psidec,
  ;           EPSIRA=epsira, EPSIDEC=epsidec, /LNP_ONLY, /NO_XYZ ])
  ;           
  ; INPUTS:
  ;       ra - Right ascension (decimal degrees). A N-elements array can be specified to calculate the Bayesian probability of several
  ;            stars at once, but then all mandatory inputs must also be N-elements arrays.
  ;       dec - Declination (decimal degrees).
  ;       pmra - Proper motion in the right ascension direction (mas/yr, must include the cos(dec) factor).
  ;       pmdec - Proper motion in the declination direction (mas/yr).
  ;       epmra - Measurement error on the proper motion in the right ascension direction (mas/yr, must not include the cos(dec) factor).
  ;       epmec -  Measurement error on the proper motion in the declination direction (mas/yr).
  ;       
  ; OPTIONAL INPUT:
  ;       PRECISION_MATRIX - Precision matrix describing the multivariate Gaussian model of the Bayesian hypothesis.
  ;                          This 6x6-elements matrix is the matrix inverse of the 6x6-elements covariance matrix, for which
  ;                          the first three dimensions describe (X,Y,Z) in parsec and the last 3 elements describe (U,V,W) in km/s.
  ;                          This keyword must be specified unless it is included within the "ASSOCIATION_STRUCTURE" keyword.
  ;                          A single PRECISION_MATRIX must be given even if N stars are analyzed at once.
  ;       CENTER_VECTOR - Central XYZUVW position of the multivariate Gaussian model (in this order; mixed units of pc and km/s).
  ;                       A single CENTER_VECTOR must be given even if N stars are analyzed at once.
  ;       PRECISION_DETERM - Determinant of the precision matrix. Specifying this keyword avoids a calculation of the determinant
  ;                          every time banyan_sigma_solve_multivar is called.
  ;       RV_MES - Radial velocity measurement to be included in the Bayesian probability (km/s).
  ;                If this keyword is set, ERV_MES must also be set.
  ;                A N-elements array must be used if N stars are analyzed at once.
  ;       ERV_MES - Measurement error on the radial velocity to be included in the Bayesian probability (km/s).
  ;                A N-elements array must be used if N stars are analyzed at once.
  ;       DIST_MES - Distance measurement to be included in the Bayesian probability (pc).
  ;                  By default, the BANYAN_SIGMA Bayesian priors are meant for this keyword to be used with trigonometric distances only.
  ;                  Otherwise, the rate of true positives may be far from the nominal values described in Gagné et al. (ApJS, 2017, XX, XX).
  ;                  If this keyword is set, EDIST_MES must also be set.
  ;                  A N-elements array must be used if N stars are analyzed at once.
  ;       EDIST_MES - Measurement error on the distance to be included in the Bayesian probability (pc).
  ;                  A N-elements array must be used if N stars are analyzed at once.
  ;       ASSOCIATION_STRUCTURE - An IDL structure that contains the keywords PRECISION_MATRIX, CENTER_VECTOR and optionally
  ;                               PRECISION_DETERM. This structure replaces the aforementioned keywords.
  ;       PSIRA - Parallax motion factor PSIRA described in Gagné et al. (ApJS, 2017, XX, XX), in units of yr^(-1).
  ;                If this keyword is set, the corresponding PSIDEC, EPSIRA and EPSIDEC keywords must also be set.
  ;                This measurement is only useful when proper motions are estimated from two single-epoch astrometric
  ;                measurements. It captures the dependence of parallax motion as a function of distance, and allows
  ;                BANYAN_SIGMA to shift the UVW center of the moving group models, which is equivalent to
  ;                correctly treating the input "proper motion" PMRA, PMDEC, EPMRA, EPMDEC as a true apparent motion.
  ;                This keyword should *not* be used if proper motions were derived from more than two epochs, or if
  ;                they were obtained from a full parallax solution. 
  ;                A N-elements array must be used if N stars are analyzed at once.
  ;       PSIDEC - Parallax motion factor PSIDEC described in Gagné et al. (ApJS, 2017, XX, XX), in units of yr^(-1).
  ;                 A N-elements array must be used if N stars are analyzed at once.
  ;       EPSIRA - Measurement error on the parallax motion factor PSIRA described in Gagné et al. (ApJS, 2017, XX, XX),
  ;                 in units of yr^(-1). A N-elements array must be used if N stars are analyzed at once.
  ;       EPSIDEC - Measurement error on the parallax motion factor PSIDEC described in Gagné et al. (ApJS, 2017, XX, XX), 
  ;                  in units of yr^(-1). A N-elements array must be used if N stars are analyzed at once.
  ;
  ; OPTIONAL INPUT KEYWORD:
  ;       /LNP_ONLY - If this keyword is set, only Bayesian probabilities will be calculated and returned.
  ;       /NO_XYZ - If this keyword is set, the width of the spatial components of the multivariate Gaussian will be widened by a large
  ;                 factor, so that the XYZ components are effectively ignored. This keyword must be used with extreme caution as it will
  ;                 generate a significant number of false-positives and confusion between the young associations.
  ;
  ; OUTPUTS:
  ;       This routine outputs an array of N IDL structures (where N is the number of stars analyzed), each of which contains:
  ;       
  ;       LN_P - Natural logarithm of the Bayesian membership probability in arbitrary units.
  ;              This probability must be normalized to mean something.
  ;       D_OPT - Optimal distance (pc) that maximizes the Bayesian likelihood for this hypothesis.
  ;       RV_OPT - Optimal radial velocity (km/s) that maximizes the Bayesian likelihood for this hypothesis.
  ;       ED_OPT - Error on the optimal distance (pc), which approximates the 68% width of how the likelihood varies with distance.
  ;       ERV_OPT - Error on the optimal radial velocity (km/s), which approximates the 68% width of how the likelihood varies with
  ;                 radial velocity.
  ;       XYZUVW - 6-dimensional array containing the XYZ and UVW position of the star at the measured radial velocity and/or
  ;                distance, or the optimal radial velocity and/or distance when the first are not available (units of pc and km/s).
  ;       EXYZUVW - Errors on XYZUVW (units of pc and km/s).
  ;       XYZ_SEP - Separation between the optimal or measured XYZ position of the star and the center of the multivariate Gaussian
  ;                 model of this Bayesian hypothesis (pc).
  ;       UVW_SEP - Separation between the optimal or measured UVW position of the star and the center of the multivariate Gaussian
  ;                 model of this Bayesian hypothesis (km/s).
  ;       XYZ_SEP - N-sigma separation between the optimal or measured XYZ position of the star and the multivariate Gaussian model
  ;                 of this Bayesian hypothesis (no units).
  ;       UVW_SEP - N-sigma separation between the optimal or measured UVW position of the star and the multivariate Gaussian model
  ;                 of this Bayesian hypothesis (no units).
  ;       MAHALANOBIS - Mahalanobis distance between the optimal or measured XYZUVW position of the star and the multivariate Gaussian
  ;                     model. A Mahalanobis distance is a generalization of a 6D N-sigma distance that accounts for covariances. 
  ;              
  ; PROCEDURES USED:
  ;       matrix_multiply_square_act, matrix_vector_product_vct, matrix_vector_product, matrix_vector_error_product, inner_product_multi, 
  ;       xyz_errors, uvw_errors, matrix_inflation_multiple, inner_product_multi_variablemetric, parabolic_cylinder_function
  ;
  ; MODIFICATION HISTORY:
  ;       WRITTEN, Jonathan Gagne, October, 25 2017
  ;-
  
  forward_function matrix_multiply_square_act, matrix_vector_product_vct, matrix_vector_product, matrix_vector_error_product, inner_product_multi, $
    xyz_errors, uvw_errors, matrix_inflation_multiple, inner_product_multi_variablemetric, parabolic_cylinder_function
  
  ;Check number of objects and array size consistency
  nobj = n_elements(ra)
  if n_elements(dec) ne nobj or n_elements(pmra) ne nobj or n_elements(pmdec) ne nobj or n_elements(epmra) ne nobj or n_elements(epmdec) ne nobj then $
    message, 'There must be the same number of each input measurement (ra,dec,pmra,pmdec,epmra,epmdec) !'
  
  ;Avoid modifying the input precision matrix
  if keyword_set(precision_matrix_in) then $
    precision_matrix = precision_matrix_in
  
  ;Read data from the young association parameters structure if it was passed on 
  if keyword_set(association_structure) then begin
    precision_matrix = association_structure.precision_matrix
    center_vector = association_structure.center_vec
    if max(strpos(strlowcase(tag_names(association_structure)),'precision_determ')) ne -1L then $
      precision_determ = association_structure.precision_determ
  endif
  
  ;Verify the array dimensions for the parallax motion factors 
  if keyword_set(psira) then begin
    if n_elements(psira) ne nobj or n_elements(psidec) ne nobj or n_elements(epsira) ne nobj or n_elements(epsidec) ne nobj then $
    message, 'There must be the same number of each input measurement (ra,dec,pmra,pmdec,epmra,epmdec,psira,psidec,epsira,epsidec) !' 
  endif
  
  ;Verify dimensions of the precision matrix
  if min((size(precision_matrix))[0L:2L] eq [2L,6L,6L]) eq 0 then $
    message, 'The covariance matrix has to have 6x6-elements !
  
  ;Verify dimensions of the center vector
  if min((size(center_vector))[0L:1L] eq [1L,6L]) eq 0 then $
    message, 'The center vector has to have 6-elements !

  ;One AU/yr in km/s divided by 1000
  kappa = 0.004743717361d0
  
  ;Compute galactic positions
  glactc_errors, ra, dec, 0d0, 0d0, gl, gb, /FAST
  
  ;Define the galactic vector
  lambda_vector = [[cos(gb*!dtor)*cos(gl*!dtor)],[cos(gb*!dtor)*sin(gl*!dtor)],[sin(gb*!dtor)]]
  
  ;Rotation Matrix for the coordinate change
  TGAL_matrix = [[-0.0548755604d0, 0.4941094279d0, -0.8676661490d0], $
    [-0.8734370902d0, -0.4448296300d0, -0.1980763734d0], $
    [-0.4838350155d0,  0.7469822445d0, 0.4559837762d0]]
  
  ;Build the A-matrix
  cra = cos(ra*!dtor)
  cdec = cos(dec*!dtor)
  sra = sin(ra*!dtor)
  sdec = sin(dec*!dtor)
  A_matrix = dblarr(nobj,3L,3L)
  A_matrix[*,0L,0L] = cra*cdec
  A_matrix[*,1L,0L] = sra*cdec
  A_matrix[*,2L,0L] = sdec
  A_matrix[*,0L,1L] = -sra
  A_matrix[*,1L,1L] = cra
  A_matrix[*,0L,2L] = -cra*sdec
  A_matrix[*,1L,2L] = -sra*sdec
  A_matrix[*,2L,2L] = cdec
  
  ;Build required vectors
  B_matrix = matrix_multiply_square_act(TGAL_matrix, A_matrix)
  M_vector = matrix_vector_product_vct(B_matrix,[1d0,0d0,0d0])
  N_vector = matrix_vector_product(B_matrix,[[replicate(0d0,nobj)],[kappa*pmra],[kappa*pmdec]])
  
  ZERO_vector = make_array(nobj,3L,value=0d0,/double)
  OMEGA_vector = [[ZERO_vector],[M_vector]]
  GAMMA_vector = [[lambda_vector],[N_vector]]
  TAU_vector = make_array(nobj,value=1d0,/double)#center_vector
  
  ;Build vectors relevant to the parallax motion factors
  PHI_vector = matrix_vector_product(B_matrix,[[replicate(0d0,nobj)],[kappa*psira],[kappa*psidec]])
  EPHI_vector = matrix_vector_error_product(B_matrix,[[replicate(0d0,nobj)],[kappa*epsira],[kappa*epsidec]])
  PSI_vector = [[ZERO_vector],[PHI_vector]]
  EPSI_vector = [[ZERO_vector],[EPHI_vector]]
  TAU_vector += PSI_vector
  
  ;Case where XYZ must be ignored
  if keyword_set(NO_XYZ) then begin
    covmat = la_invert(precision_matrix)
    covmat[0:2,0:2] = diag_matrix(dblarr(3)+1d9)
    covmat[0:2,3:*] = 0d0
    covmat[3:*,0:2] = 0d0
    precision_matrix = la_invert(covmat)
  endif
  
  ;Calculate the BANYAN moments
  OMEGA_OMEGA = inner_product_multi(OMEGA_vector, OMEGA_vector, precision_matrix)
  GAMMA_GAMMA = inner_product_multi(GAMMA_vector, GAMMA_vector, precision_matrix)
  OMEGA_GAMMA = inner_product_multi(OMEGA_vector, GAMMA_vector, precision_matrix)
  OMEGA_TAU = inner_product_multi(OMEGA_vector, TAU_vector, precision_matrix)
  GAMMA_TAU = inner_product_multi(GAMMA_vector, TAU_vector, precision_matrix)
  TAU_TAU = inner_product_multi(TAU_vector, TAU_vector, precision_matrix)

  ;Include RV and distance measurements
  if keyword_set(dist_mes) then begin
    good = where(finite(dist_mes), ngood)
    if ngood ne 0L then begin
      GAMMA_GAMMA[good] += 1d0/(edist_mes[good]>1d-3)^2
      GAMMA_TAU[good] += dist_mes[good]/(edist_mes[good]>1d-3)^2
      TAU_TAU[good] += dist_mes[good]^2/(edist_mes[good]>1d-3)^2
    endif
  endif
  if keyword_set(rv_mes) then begin
    good = where(finite(rv_mes), ngood)
    if ngood ne 0L then begin
      OMEGA_OMEGA[good] += 1d0/(erv_mes[good]>1d-3)^2
      OMEGA_TAU[good] += rv_mes[good]/(erv_mes[good]>1d-3)^2
      TAU_TAU[good] += rv_mes[good]^2/(erv_mes[good]>1d-3)^2
    endif
  endif

  ;Calculate determinant of precision matrix if not already specified
  if ~keyword_set(precision_determ) then $
    determ_precision = la_determ(precision_matrix) else $
    determ_precision = precision_determ
  determ_precision0 = determ_precision
  if determ_precision le 0 then message, 'Ill-defined value for the precision matrix determinant (covariance matrix has <= 0 determinant) !'

  ;Calculate first estimates for the optimal distance and RV (this estimate ignores error bars on proper motion)
  beta = (GAMMA_GAMMA-OMEGA_GAMMA^2/OMEGA_OMEGA)/2d0
  if min(beta,/nan) lt 0 then message, 'Ill-defined value for beta !'
  gamma = (OMEGA_GAMMA*OMEGA_TAU/OMEGA_OMEGA-GAMMA_TAU)
  d_opt = (sqrt(gamma^2+32d0*beta)-gamma)/(4d0*beta)
  rv_opt = (4d0-GAMMA_GAMMA*d_opt^2+GAMMA_TAU*d_opt)/(OMEGA_GAMMA*d_opt)

  ;Propagate measurement errors if they are set
  determ_precision = replicate(determ_precision,nobj)
  if determ_precision[0L] lt 0 then $
    message, 'Ill-defined value for the precision matrix determinant (covariance matrix has <= 0 determinant) !'
  
  ;Propagate measurement errors to XYZ and UVW
  xyz_errors, ra, dec, d_opt, 0d0, 0d0, 0d0, X, Y, Z, EX, EY, EZ
  uvw_errors, ra, dec, pmra, pmdec, rv_opt, d_opt, 0d0, 0d0, epmra, epmdec, 0d0, 0d0, U, V, W, EU, EV, EW
  
  ;Define a vector of inflation for the covariance matrix
  inflate_vec = [[EX],[EY],[EZ],[EU],[EV],[EW]]
  
  ;Include the error bars on PSI_VEC. Since it shifts the moving group UVW centers, its error bars are added in quadrature here
  inflate_vec = sqrt(inflate_vec^2d0 + PSI_vector^2d0)
  
  covariance_matrix = la_invert(precision_matrix)
  diag_covariance = diag_matrix(covariance_matrix)
  
  ;Apply diagonal inflation on the covariance matrix while avoiding a for loop
  ones_vector1 = make_array(nobj,value=1d0,/double)
  inflation_factors = 1d0 + inflate_vec^2/(ones_vector1#diag_covariance)
  inflation_determ = exp(total(alog(inflation_factors),2L,/nan))
  
  ;Avoid singular matrices
  bad = where(inflation_determ le 0, nbad)
  if nbad ne 0 then message, 'Ill-defined value for the precision matrix determinant (covariance matrix has <= 0 determinant) !'
  
  ;Inflate the precision matrices
  precision_matrix_inflated = matrix_inflation_multiple(precision_matrix, 1d0/sqrt(inflation_factors))
  
  ;Calculate new determinants for the now multiple precision matrices
  determ_precision = determ_precision0/inflation_determ

  ;Update the BANYAN moments
  OMEGA_OMEGA = inner_product_multi_variablemetric(OMEGA_vector, OMEGA_vector, precision_matrix_inflated)
  GAMMA_GAMMA = inner_product_multi_variablemetric(GAMMA_vector, GAMMA_vector, precision_matrix_inflated)
  OMEGA_GAMMA = inner_product_multi_variablemetric(OMEGA_vector, GAMMA_vector, precision_matrix_inflated)
  OMEGA_TAU = inner_product_multi_variablemetric(OMEGA_vector, TAU_vector, precision_matrix_inflated)
  GAMMA_TAU = inner_product_multi_variablemetric(GAMMA_vector, TAU_vector, precision_matrix_inflated)
  TAU_TAU = inner_product_multi_variablemetric(TAU_vector, TAU_vector, precision_matrix_inflated)
  
  ;Include RV and distance measurements and uncertainties, where applicable
  if keyword_set(dist_mes) then begin
    good = where(finite(dist_mes), ngood)
    if ngood ne 0L then begin
      GAMMA_GAMMA[good] += 1d0/edist_mes[good]^2
      GAMMA_TAU[good] += dist_mes[good]/edist_mes[good]^2
      TAU_TAU[good] += dist_mes[good]^2/edist_mes[good]^2
    endif
  endif
  if keyword_set(rv_mes) then begin
    good = where(finite(rv_mes), ngood)
    if ngood ne 0L then begin
      OMEGA_OMEGA[good] += 1d0/erv_mes[good]^2
      OMEGA_TAU[good] += rv_mes[good]/erv_mes[good]^2
      TAU_TAU[good] += rv_mes[good]^2/erv_mes[good]^2
    endif
  endif

  ;Calculate more relevant quantities
  beta = (GAMMA_GAMMA-OMEGA_GAMMA^2/OMEGA_OMEGA)/2d0
  gamma = (OMEGA_GAMMA*OMEGA_TAU/OMEGA_OMEGA-GAMMA_TAU)
  if min(beta,/nan) lt 0 then message, 'Ill-defined value for beta !'
  
  ;Update optimal distance and RV
  d_opt = (sqrt(gamma^2d0+32d0*beta)-gamma)/(4d0*beta)
  rv_opt = (4d0-GAMMA_GAMMA*d_opt^2+GAMMA_TAU*d_opt)/(OMEGA_GAMMA*d_opt)
  erv_opt = 1d0/sqrt(OMEGA_OMEGA)
  ed_opt = 1d0/sqrt(GAMMA_GAMMA)
  
  ;Calculate probabilities
  zeta = (TAU_TAU-OMEGA_TAU^2/OMEGA_OMEGA)/2d0
  xarg = gamma/sqrt(2d0*beta)
  ln_P_coeff = -0.5d0*alog(OMEGA_OMEGA)-2.5d0*alog(beta)+0.5d0*alog(determ_precision)
  ln_P_part1 = xarg^2/2d0 - zeta
  ln_P_part2 = alog(parabolic_cylinder_function(-5L, xarg, /FPRIME)>1d-318)
  ln_P = ln_P_coeff + ln_P_part1 + ln_P_part2
  
  ;Return the probability if this is the only required output
  if keyword_set(lnp_only) then return, ln_P
  
  ;Use the true distance and/or radial velocity if they were specified
  dist_final = d_opt
  edist_final = ed_opt
  if keyword_set(dist_mes) then begin
    gdist = where(finite(dist_mes), ngdist)
    if ngdist ne 0L then begin
      dist_final[gdist] = dist_mes[gdist]
      edist_final[gdist] = edist_mes[gdist]
    endif
  endif
  rv_final = rv_opt
  erv_final = erv_opt
  if keyword_set(rv_mes) then begin
    grv = where(finite(rv_mes), ngrv)
    if ngrv ne 0L then begin
      rv_final[grv] = rv_mes[grv]
      erv_final[grv] = erv_mes[grv]
    endif
  endif
  
  ;Calculate optimal XYZUVW positions while ignoring error bars on ra and dec 
  xyz_errors, ra, dec, dist_final, finite(ra)*0., finite(dec)*0., edist_final, XO, YO, ZO, EXO, EYO, EZO
  uvw_errors, ra, dec, pmra, pmdec, rv_final, $
    dist_final, finite(ra)*0., finite(dec)*0., epmra, epmdec, $
    erv_final, edist_final, UO, VO, WO, EUO, EVO, EWO
  XYZUVW = [[XO],[YO],[ZO],[UO],[VO],[WO]]
  EXYZUVW = [[EXO],[EYO],[EZO],[EUO],[EVO],[EWO]]
  
  ;Calculate mahalanobis distance between optimal position and the multivariate Gaussian model
  vec = XYZUVW - TAU_vector
  mahalanobis = sqrt(inner_product_multi_variablemetric(vec, vec, precision_matrix_inflated))
  
  XYZ_sep = sqrt(total((XYZUVW[*,0:2]-TAU_vector[*,0:2])^2,2,/nan))
  UVW_sep = sqrt(total((XYZUVW[*,3:*]-TAU_vector[*,3:*])^2,2,/nan))
  
  XYZ_sig = sqrt(inner_product_multi_variablemetric(vec[*,0:2], vec[*,0:2], precision_matrix_inflated[*,0:2,0:2]))
  UVW_sig = sqrt(inner_product_multi_variablemetric(vec[*,3:5], vec[*,3:5], precision_matrix_inflated[*,3:5,3:5]))
  
  ;Create an output structure
  output_structure = {LN_P:!values.d_nan, D_OPT:!values.d_nan, RV_OPT:!values.d_nan, ED_OPT:!values.d_nan, ERV_OPT:!values.d_nan, $
    XYZUVW:dblarr(6)+!values.d_nan, EXYZUVW:dblarr(6)+!values.d_nan, XYZ_SEP:!values.d_nan, UVW_SEP:!values.d_nan, $
    XYZ_SIG:!values.d_nan, UVW_SIG:!values.d_nan, MAHALANOBIS:!values.d_nan}
  output_structure = replicate(output_structure, nobj)
  
  ;Fill up the output structure
  output_structure.LN_P = ln_P
  output_structure.D_OPT = d_opt
  output_structure.RV_OPT = rv_opt
  output_structure.ED_OPT = ed_opt
  output_structure.ERV_OPT = erv_opt
  output_structure.XYZUVW = transpose(xyzuvw)
  output_structure.EXYZUVW = transpose(exyzuvw)
  output_structure.XYZ_SEP = xyz_sep
  output_structure.UVW_SEP = uvw_sep
  output_structure.XYZ_SIG = xyz_sig
  output_structure.UVW_SIG = uvw_sig
  output_structure.MAHALANOBIS = mahalanobis
  
  ;Return the output structure
  return, output_structure
  
End