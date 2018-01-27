Pro uvw_errors, ra, dec, pmra, pmdec, vrad, dist, e_ra, e_dec, e_pmra, e_pmdec, e_vrad, e_dist, $
                U, V, W, E_U, E_V, E_W, PLX=plx, FAST=fast, PRECISE=precise, NOERRORS=noerrors, $ 
                TRANSFORMATION_MATRIX=transformation_matrix, TRANSFORMATION_VECTOR=transformation_vector, $
                E_TRANSFORMATION_VECTOR=e_transformation_vector
  ;+
  ; NAME:
  ;       UVW_ERRORS
  ;
  ; PURPOSE:
  ;       Computes space velocities (U,V,W) and their errors (EU,EV,EW) from astrometry (ra, dec), distance (in parsec), proper motion (pmra, pmdec in mas/yr), radial velocity (in km/s) and their errors (same units).
  ;       U points towards the galactic center, V towards the local direction of rotation in the plane of the galaxy, and W towards the North Galactic Pole
  ;       (so that the UVW system is right-handed).
  ;
  ; CALLING SEQUENCE:
  ;       UVW_ERRORS, ra, dec, pmra, pmdec, vrad, dist, e_ra, e_dec, e_pmra, e_pmdec, e_vrad, e_dist, $
  ;                U, V, W, E_U, E_V, E_W[, PLX=plx, LSR=lsr, FAST=fast, PRECISE=precise]
  ;
  ; INPUTS:
  ;       RA = Scalar or 1D array indicating right ascension (in decimal degrees).
  ;       DEC = Scalar or 1D array indicating declination (in decimal degrees).
  ;       PMRA = Scalar or 1D array indicating proper motion in the direction of right ascension (in mas/yr).
  ;       PMDEC = Scalar or 1D array indicating proper motion in the direction of declination (in mas/yr).
  ;       VRAD = Scalar or 1D array indicating radial velocity (in km/s).
  ;       DIST = Scalar or 1D array indicating distance (in parsec).
  ;       ERA = Scalar or 1D array indicating the error on the right ascension (in decimal degrees) - These are not considered unless the /PRECISE keyword is used.
  ;       EDEC = Scalar or 1D array indicating the error on the declination (in decimal degrees) - These are not considered unless the /PRECISE keyword is used.
  ;       EPMRA = Scalar or 1D array indicating the error on proper motion in the direction of right ascension (in mas/yr).
  ;       EPMDEC = Scalar or 1D array indicating the error on proper motion in the direction of declination (in mas/yr).
  ;       VRAD = Scalar or 1D array indicating the error on radial velocity (in km/s).
  ;       EDIST = Scalar or 1D array indicating the error on distance (in parsec).
  
  ; OPTIONAL INPUTS:
  ;       PLX = Scalar or 1D array indicating parallax (in mas ; overrides the DIST input).
  ;       EPLX = Scalar or 1D array indicating the error on parallax (in mas ; overrides the DIST input).
  ;
  ; OPTIONAL INPUT KEYWORD:
  ;       /FAST - If this keyword is set, then no error propagation is computed.
  ;       /PRECISE - If this keyword is set, then errors on RA and DEC are propagated, and a statistical nonlinear analysis (including covariance factors) is performed.
  ;
  ; OUTPUTS:
  ;       X, Y, Z = Heliocentric galactic position. See "Purpose" for more details".
  ;       EX, EY, EZ = Errors on heliocentric galactic position. See "Purpose" for more details".
  ;
  ; NOTES :
  ;       This routine is based on equations from gal_uvw.pro (astrolib), except for the analytical propagation of errors.
  ;
  ; RESTRICTIONS:
  ;       (1) If you give 1D arrays as inputs, each of them must have the same size.
  ;
  ; PROCEDURES USED:
  ;       NONE
  ;
  ; MODIFICATION HISTORY:
  ;       WRITTEN, Jonathan Gagne, October 10th, 2012.
  ;-
  
  ;Transform parallax into distances when specified
  if keyword_set(plx) then begin
    if ~keyword_set(fast) then $
      e_plx = e_dist
    plx = dist
  endif else begin
    if ~keyword_set(fast) then $
      e_plx = e_dist*1d3/dist^2
    plx = 1d3/dist
  endelse
  
  ;Rotation Matrix for the coordinate change
  ;The first column is inverted with respect to the Astrolib convention, such that
  ;U points toward the Galactic center and UVW forms a right-handed system.
  A = [[-0.0548755604d0, +0.4941094279d0, -0.8676661490d0], $
       [-0.8734370902d0, -0.4448296300d0, -0.1980763734d0], $
       [-0.4838350155d0,  0.7469822445d0, +0.4559837762d0]]
   
  ;1 A.U/yr in km/s
  k = 4.743717361d0
  
  ;Shortcuts
  cosd = cos(dec*!dtor)
  sind = sin(dec*!dtor)
  cosa = cos(ra*!dtor)
  sina = sin(ra*!dtor)
  T1 =  A[0,0]*cosa*cosd + A[0,1]*sina*cosd + A[0,2]*sind
  T2 = -A[0,0]*sina      + A[0,1]*cosa
  T3 = -A[0,0]*cosa*sind - A[0,1]*sina*sind + A[0,2]*cosd
  T4 =  A[1,0]*cosa*cosd + A[1,1]*sina*cosd + A[1,2]*sind
  T5 = -A[1,0]*sina      + A[1,1]*cosa
  T6 = -A[1,0]*cosa*sind - A[1,1]*sina*sind + A[1,2]*cosd
  T7 =  A[2,0]*cosa*cosd + A[2,1]*sina*cosd + A[2,2]*sind
  T8 = -A[2,0]*sina      + A[2,1]*cosa                   
  T9 = -A[2,0]*cosa*sind - A[2,1]*sina*sind + A[2,2]*cosd
  
  ;Compute analytical derivatives to propagate error bars
  if ~keyword_set(fast) and ~keyword_set(precise) then begin
    T1_alpha =  T2*cosd
    T2_alpha = -A[0,0]*cosa      - A[0,1]*sina
    T3_alpha =  A[0,0]*sina*sind - A[0,1]*cosa*sind
    T4_alpha =  T5*cosd
    T5_alpha = -A[1,0]*cosa      - A[1,1]*sina
    T6_alpha =  A[1,0]*sina*sind - A[1,1]*cosa*sind
    T7_alpha =  T8*cosd
    T8_alpha = -A[2,0]*cosa      - A[2,1]*sina
    T9_alpha =  A[2,0]*sina*sind - A[2,1]*cosa*sind
    T1_delta =  T3
    T2_delta =  0d0
    T3_delta = -T1
    T4_delta =  T6
    T5_delta =  0d0
    T6_delta = -T4
    T7_delta =  T9
    T8_delta =  0d0
    T9_delta = -T7
  endif
  
  ;Compute U,V,W
  U = T1*vrad + T2*k*pmra/plx + T3*k*pmdec/plx
  V = T4*vrad + T5*k*pmra/plx + T6*k*pmdec/plx
  W = T7*vrad + T8*k*pmra/plx + T9*k*pmdec/plx
  
  ;If only one star is analyzed, also calculate transormation vector and matrix
  ; which are useful to display error bars on UVW correctly
  if n_elements(T1) eq 1L then begin
    transformation_matrix = [[t1,t2,t3],[t4,t5,t6],[t7,t8,t9]]
    transformation_vector = [vrad, k*pmra/plx, k*pmdec/plx]
  endif
  
  ;Return now if no errurs are specified
  if keyword_set(noerrors) then return
  
  ;Propagate errors with analytical approximations
  if ~keyword_set(fast) then begin
    iplx = 1d0/plx
    e_iplx = 1d0/plx^2*e_plx
    if n_elements(t1) eq 1L then $
      e_transformation_vector = [e_vrad, k*sqrt((e_pmra/plx)^2+(pmra*e_iplx)^2), k*sqrt((e_pmdec/plx)^2+(pmdec*e_iplx)^2)]
    if keyword_set(precise) then begin
      E_U = sqrt( T1^2*e_vrad^2 + T2^2*k^2* ( e_pmra^2 * iplx^2 + e_iplx^2 * pmra^2 + e_pmra^2 * e_iplx^2 ) + $  
                  T3^2*k^2* ( e_pmdec^2 * iplx^2 + e_iplx^2 * pmdec^2 + e_pmdec^2 * e_iplx^2 ) )
      E_V = sqrt( T4^2*e_vrad^2 + T5^2*k^2* ( e_pmra^2 * iplx^2 + e_iplx^2 * pmra^2 + e_pmra^2 * e_iplx^2 ) + $  
                  T6^2*k^2* ( e_pmdec^2 * iplx^2 + e_iplx^2 * pmdec^2 + e_pmdec^2 * e_iplx^2 ) )
      E_W = sqrt( T7^2*e_vrad^2 + T8^2*k^2* ( e_pmra^2 * iplx^2 + e_iplx^2 * pmra^2 + e_pmra^2 * e_iplx^2 ) + $  
                  T9^2*k^2* ( e_pmdec^2 * iplx^2 + e_iplx^2 * pmdec^2 + e_pmdec^2 * e_iplx^2 ) )
    endif else begin
      E_U = sqrt( (T1*e_vrad)^2 + (k*T2/plx*e_pmra)^2 + (k*T3/plx*e_pmdec)^2 + (T1_alpha*vrad + T2_alpha*k*pmra/plx + $
                   T3_alpha*k*pmdec/plx)^2*e_ra^2 + (T1_delta*vrad + T2_delta*k*pmra/plx + T3_delta*k*pmdec/plx)^2*e_dec^2 + $
                  (k/plx^2*(pmra*T2+pmdec*T3)*e_plx)^2 )
      E_V = sqrt( (T4*e_vrad)^2 + (k*T5/plx*e_pmra)^2 + (k*T6/plx*e_pmdec)^2 + (T4_alpha*vrad + T5_alpha*k*pmra/plx + $
                   T6_alpha*k*pmdec/plx)^2*e_ra^2 + (T4_delta*vrad + T5_delta*k*pmra/plx + T6_delta*k*pmdec/plx)^2*e_dec^2 + $
                  (k/plx^2*(pmra*T5+pmdec*T6)*e_plx)^2 )
      E_W = sqrt( (T7*e_vrad)^2 + (k*T8/plx*e_pmra)^2 + (k*T9/plx*e_pmdec)^2 + (T7_alpha*vrad + T8_alpha*k*pmra/plx + $
                   T9_alpha*k*pmdec/plx)^2*e_ra^2 + (T7_delta*vrad + T8_delta*k*pmra/plx + T9_delta*k*pmdec/plx)^2*e_dec^2 + $
                  (k/plx^2*(pmra*T8+pmdec*T9)*e_plx)^2 )
    endelse
  endif else begin
    e_u = finite(U)*0.
    e_v = e_u
    e_w = e_u
  endelse
  
End