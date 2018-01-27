Function uniq_unsorted, vec, INDFULL=indfull, bad=bad
  ;Written by J. Gagne to return unique indices without sorting data
  
  forward_function remove
  
  uvec = vec
  uvec = uvec[sort(uvec)]
  uvec = uvec[uniq(uvec)]
  nu = n_elements(uvec)
  nvec = n_elements(vec)
  uvec_used = intarr(nu)
  
  ind = intarr(nvec)-1
  for i=0L, nu-1L do begin
    gg = where(vec eq uvec[i], ngg)
    if ngg eq 0L then continue
    ind[gg[0L]] = gg[0]
  endfor
  indfull = ind
  bad = where(ind eq -1L, nbad)
  if nbad ne 0L then remove, bad, ind
  return, ind
End