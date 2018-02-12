Pro xyz_errors, ra, dec, dist, era, edec, edist, X, Y, Z, EX, EY, EZ, FAST=fast
;+
; NAME:
;       XYZ_ERRORS
;       
; PURPOSE:
;       Computes heliocentric galactic coordinates (X,Y,Z) and their errors (EX,EY,EZ) from astrometry (ra, dec) and distance (in parsec) and their errors.
;       X points towards the galactic center, Y towards the local direction of rotation in the plane of the galaxy, and Z towards the North Galactic Pole
;       (so that the XYZ system is right-handed).
;       
; CALLING SEQUENCE:
;       XYZ_ERRORS, ra, dec, dist, era, edec, edist, X, Y, Z[, EX, EY, EZ, /FAST]
;
; INPUTS:
;       RA = Scalar or 1D array indicating right ascension (in decimal degrees).
;       DEC = Scalar or 1D array indicating declination (in decimal degrees).
;       DIST = Scalar or 1D array indicating distance (in parsec).
;       ERA = Scalar or 1D array indicating the error on the right ascension (in decimal degrees).
;       EDEC = Scalar or 1D array indicating the error on the declination (in decimal degrees).
;       EDIST = Scalar or 1D array indicating the error on distance (in parsec).
;
; OPTIONAL INPUT KEYWORD:
;       /FAST - If this keyword is set, then no error propagation is computed.
;
; OUTPUTS:
;       X, Y, Z = Heliocentric galactic position. See "Purpose" for more details".
;       EX, EY, EZ = Errors on heliocentric galactic position. See "Purpose" for more details".
; NOTES :
;       This routine is based on equations from gal_uvw.pro (astrolib), except for the analytical propagation of errors.
;
; RESTRICTIONS:
;       (1) If you give 1D arrays as inputs, each of them must have the same size.
;
; PROCEDURES USED:
;       GLACTC_ERRORS
;
; MODIFICATION HISTORY:
;       WRITTEN, Jonathan Gagne, October 10th, 2012.
;-
  
  ;Declare subroutines
  forward_function glactc_errors
  
  ;Transform (RA,DEC) into Galactic coordinates
  glactc_errors, ra, dec, era, edec, gl, gb, egl, egb, /FAST;=fast
  egl = gl*0d0
  egb = gb*0d0
  
  ;Shortcuts
  cgb = cos(gb*!dtor)
  sgb = sin(gb*!dtor)
  cgl = cos(gl*!dtor)
  sgl = sin(gl*!dtor)
  
  ;Compute Galactic positions X,Y,Z 
  X = (cgb*cgl)*dist
  Y = (cgb*sgl)*dist
  Z = sgb*dist
  
  ;Propagate errors with analytical approximations
  if ~keyword_set(fast) then begin
    EX = sqrt(((sgb*cgl*egb)^2+(cgb*sgl*egl)^2)*(!dtor*dist)^2+(cgb*cgl*edist)^2)
    EY = sqrt(((sgb*sgl*egb)^2+(cgb*cgl*egl)^2)*(!dtor*dist)^2+(cgb*sgl*edist)^2)
    EZ = sqrt(((cgb*egb)*dist*!dtor)^2+(sgb*edist)^2)
  endif else begin
    ex = finite(x)*0.
    ey = ex
    ez = ex
  endelse
  
End