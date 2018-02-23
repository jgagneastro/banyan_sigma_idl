Function banyan_sigma_cmd_to_distance, data_file, xmag1, xmag2, ymag, exmag1, exmag2, eymag
  ;This function returns the expected distance for each star, and each age category, given its
  ;  position in a absolute Gaia G mag versus Gaia-2MASS G-J color

  if ~keyword_set(data_file) then $
    message, ' A FITS file with the CMD data must be specified !'
  
  if ~file_test(data_file) then $
    message, ' The CMD data file could not be found !'
  
  ;Determine number of stars
  nstars = n_elements(xmag1)

  ;Restore CMD data
  cmd_data = mrdfits(data_file,1,/silent)

  ;Determine colors
  colors = xmag1 - xmag2

  ;Determine error bar on colors
  e_colors = sqrt(exmag1^2 + exmag2^2)

  ;Loop on classes of young associations to determine photometric distances
  nseqs = n_elements(cmd_data.NAMES)
  group_name_separator = '|'
  delta_mag_seqs = dblarr(nstars,nseqs)+!values.d_nan
  e_delta_mag_seqs = dblarr(nstars,nseqs)+!values.d_nan
  all_group_names = []
  for i=0L, nseqs-1L do begin

    ;Store group names
    all_group_names = [all_group_names, strsplit(strtrim(cmd_data.NAMES[i],2),'|',/extract)]

    ;Determine expected absolute G mag at the G-J colors of targets
    abs_mag = interpol(cmd_data.SEQS_Y[*,i], cmd_data.SEQS_X, colors)

    ;Determine deviations in expected absolute G mags given error bar on color
    e_color_contribution = interpol(deriv(cmd_data.SEQS_X,cmd_data.SEQS_Y[*,i]), cmd_data.SEQS_X, colors)*e_colors

    ;Determine error bar
    e_abs_mag_neg = interpol(reform(cmd_data.SEQS_ERR[*,0,i]), cmd_data.SEQS_X, colors)
    e_abs_mag_pos = interpol(reform(cmd_data.SEQS_ERR[*,1,i]), cmd_data.SEQS_X, colors)
    e_abs_mag = (e_abs_mag_neg + e_abs_mag_pos)/2d0

    ;Add measured error bars on gmag and color in quadrature
    e_abs_mag = sqrt(e_abs_mag^2 + e_color_contribution^2)

    delta_mag_seqs[*,i] = ymag - abs_mag
    e_delta_mag_seqs[*,i] = sqrt(eymag^2 + e_abs_mag^2)

  endfor

  ;Determine photometric distances
  phot_dist = 10d0^(delta_mag_seqs/5d0+1)
  e_phot_dist = alog(10d0)*phot_dist/5d0*e_delta_mag_seqs

  ;Remove distance estimates for objects with G-J colors outside the range
  bad = where(colors lt min(cmd_data.SEQS_X,/nan) or colors gt max(cmd_data.SEQS_X,/nan), nbad)
  if nbad ne 0L then begin
    phot_dist[bad,*] = !values.d_nan
    e_phot_dist[bad,*] = !values.d_nan
  endif

  ;Create an output structure for BANYAN Sigma
  ngroups = n_elements(all_group_names)
  create_struct, out, '', all_group_names, strjoin(strarr(ngroups)+'D',',')

  ;Replicate structure
  nan_str, out
  out = replicate(out, nstars)
  out_error = out

  ;Fill the structure
  out_tags = tag_names(out)
  for i=0L, ngroups-1L do begin
    gbatch = where(strpos(cmd_data.NAMES, group_name_separator+all_group_names[i]+group_name_separator) ne -1L, ngbatch)
    if ngbatch ne 1 then stop
    gtag = where(out_tags eq all_group_names[i], ngtag)
    if ngtag ne 1 then stop

    out.(gtag[0L]) = phot_dist[*,gbatch[0L]]
    out_error.(gtag[0L]) = e_phot_dist[*,gbatch[0L]]
  endfor
  
;  bad = where(~finite(phot_dist) and finite(e_phot_dist), nbad)
;  if nbad ne 0L then stop
;  bad = where(finite(phot_dist) and ~finite(e_phot_dist), nbad)
;  if nbad ne 0L then stop

  return, {CONSTRAINT_DIST_PER_HYP:out,CONSTRAINT_EDIST_PER_HYP:out_error}

End