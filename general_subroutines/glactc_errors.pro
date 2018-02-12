Pro glactc_errors, ra, dec, era, edec, gl, gb, egl, egb, FAST=fast
;+
; NAME:
;       GLACTC_ERRORS
;
; PURPOSE:
;       This procedure is similar to glactc.pro (j=1 direction) from the astrolib, but it takes measurement errors into account. Also, it does not need to
;       precess coordinates (which significantly fastens cases with many objects).
;
; CALLING SEQUENCE:
;       GLACTC_ERRORS, ra, dec, era, edec, gl, gb[, egl, egb, /FAST]
;
; INPUTS:
;       RA = Scalar or 1D array indicating right ascension (in decimal degrees).
;       DEC = Scalar or 1D array indicating declination (in decimal degrees).
;       ERA = Scalar or 1D array indicating the error on the right ascension (in decimal degrees).
;       EDEC = Scalar or 1D array indicating the error on the declination (in decimal degrees).
;
; OPTIONAL INPUT KEYWORD:
;       /FAST - If this keyword is set, then no error propagation is computed.
;
; OUTPUTS:
;       GL, GB = Galactic longitude and latitude, respectively (in degrees).
;       EGL, EGB = Errors on Galactic longitude and latitude, respectively (in degrees).

; NOTES :
;       This routine is based on equations from glactc.pro (astrolib), except they were 
;       modified to avoid the need of precession, and analytical propagation of errors
;       were carried and added by the authors.
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
  
  radeg = 1.8d2/!dpi
  
  ;Avoid calculating error bars with /FAST
  if keyword_set(fast) then begin
    gl = double(!values.f_nan)
    gb = gl
  endif
  
  ;J2000.0 positions of the Galactic North pole (b=90 degrees), from Carrol and Ostlie.
  rapol = 192.8595d0 
  decpol = 27.12825d0
  
  ;Latitude of the Celestial Galactic North Pole (delta = 90 degrees at J2000.0).
  ;There is an error in Carrol & Ostlie for this value (they list 123d55m55.2s at p.900), which is inconsistent with the
  ;NASA astrolib value in glactc.pro (122d55m55.2s), which we adopt here. 
  lnorth = 122.932d0
  
  ;Shortcuts
  ra2 = ra - rapol
  sdp = sin(decpol/radeg)
  cdp = cos(decpol/radeg)
  sra = sin(ra2/radeg)
  cra = cos(ra2/radeg)
  sdec = sin(dec/radeg)
  cdec = cos(dec/radeg)
  
  ;Calculate Galactic latitude
  gamma = sdp*sdec + cdp*cdec*cra
  gb = radeg*asin(gamma)
  ;Analytical derivatives from Wolfram Alpha are used to propagate error bars
  if ~keyword_set(fast) then begin
    der_den = sqrt(1-(cdp*cra*cdec+sdp*sdec)^2)
    gb_dec = (sdp*cdec - cdp*cra*sdec)/der_den
    gb_ra = cdp*sra*cdec/der_den
    egb = sqrt((gb_dec*edec)^2+(gb_ra*era)^2)
  endif
  
  ;Calculate Galactic longitude
  gl = lnorth - radeg*atan((cdec*sra),(sdec-sdp*gamma)/cdp)
  gl = (gl + 360d0) mod 360d0
  ;Analytical derivatives from Wolfram Alpha are used to propagate error bars
  if ~keyword_set(fast) then begin
    der_den = (cdec^2*(sdp^2*cra^2+sra^2)-2*sdp*cdp*cra*sdec*cdec+cdp^2*sdec^2)
    gl_ra = cdec*(cdp*cra*sdec-sdp*cdec)/der_den
    gl_dec = cdp*sra/der_den
    egl = sqrt((gl_ra*era)^2+(gl_dec*edec)^2)
  endif
  
End