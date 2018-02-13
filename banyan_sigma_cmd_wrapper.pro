Function banyan_sigma_cmd_wrapper, stars_data, DATA_FILE=data_file, column_names=column_names, _EXTRA=extra, $
  CONSTRAINT_DIST_PER_HYP=constraint_dist_per_hyp, CONSTRAINT_EDIST_PER_HYP=constraint_edist_per_hyp
  ;COLUMN_NAMES must contain XMAG1, XMAG2 and YMAG, EXMAG1, EXMAG2 and EYMAG
  
  ;Default data file
  if ~keyword_set(data_file) then $
    data_file = 'cmd_tgas_gj_g.fits'
  
  ;Add directory if needed and verify if data file exists
  if file_dirname(data_file) eq '.' then $
    data_file = file_dirname(routine_filepath())+path_sep()+'data'+path_sep()+data_file
  if ~file_test(data_file) then $
    message, ' The CMD data file could not be found !'
  
  ;Verify tags of column_names
  column_tags = tag_names(column_names)
  if ~max(column_tags eq 'XMAG1') then $
    message, ' The tag XMAG1 must be specified in COLUMN_NAMES !'
  if ~max(column_tags eq 'EXMAG1') then $
    message, ' The tag EXMAG1 must be specified in COLUMN_NAMES !'
  if ~max(column_tags eq 'XMAG2') then $
    message, ' The tag XMAG2 must be specified in COLUMN_NAMES !'
  if ~max(column_tags eq 'EXMAG2') then $
    message, ' The tag EXMAG2 must be specified in COLUMN_NAMES !'
  if ~max(column_tags eq 'YMAG') then $
    message, ' The tag YMAG must be specified in COLUMN_NAMES !'
  if ~max(column_tags eq 'EYMAG') then $
    message, ' The tag EYMAG must be specified in COLUMN_NAMES !'
  
  ;Find the correct tags in the stars data
  stars_tags = tag_names(stars_data)
  xmag1tag = where(stars_tags eq column_names.XMAG1, ntag)
  if ntag ne 1 then $
    message, ' The tag '+strtrim(column_names.XMAG1)+' was not specified in the input structure !'
  exmag1tag = where(stars_tags eq column_names.EXMAG1, ntag)
  if ntag ne 1 then $
    message, ' The tag '+strtrim(column_names.EXMAG1)+' was not specified in the input structure !'
  xmag2tag = where(stars_tags eq column_names.XMAG2, ntag)
  if ntag ne 1 then $
    message, ' The tag '+strtrim(column_names.XMAG2)+' was not specified in the input structure !'
  exmag2tag = where(stars_tags eq column_names.EXMAG2, ntag)
  if ntag ne 1 then $
    message, ' The tag '+strtrim(column_names.EXMAG2)+' was not specified in the input structure !'
  ymagtag = where(stars_tags eq column_names.YMAG, ntag)
  if ntag ne 1 then $
    message, ' The tag '+strtrim(column_names.YMAG)+' was not specified in the input structure !'
  eymagtag = where(stars_tags eq column_names.EYMAG, ntag)
  if ntag ne 1 then $
    message, ' The tag '+strtrim(column_names.EYMAG)+' was not specified in the input structure !'
  
  ;Call the CMD => Distance function
  cmd_output = banyan_sigma_cmd_to_distance(data_file, stars_data.(xmag1tag[0L]), stars_data.(xmag2tag[0L]), stars_data.(ymagtag[0L]), stars_data.(exmag1tag[0L]), stars_data.(exmag2tag[0L]), stars_data.(eymagtag[0L]))
  
  constraint_dist_per_hyp = cmd_output.CONSTRAINT_DIST_PER_HYP
  constraint_edist_per_hyp = cmd_output.CONSTRAINT_EDIST_PER_HYP
  
  ;Call BANYAN Sigma and provide it with the CMD data
  outs = banyan_sigma(stars_data, COLUMN_NAMES=column_names, CONSTRAINT_DIST_PER_HYP=constraint_dist_per_hyp, $
    CONSTRAINT_EDIST_PER_HYP=constraint_edist_per_hyp, _extra=extra)
  
  ;Return outputs to caller
  return, outs
  
End