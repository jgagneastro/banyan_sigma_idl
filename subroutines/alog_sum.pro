Function alog_sum, x
  ;This function was written by J. Gagne to sum a 1D array X along a given axis in log space.
  ;X must already be a natural log quantity.
  
  ;Find a normalization factor that avoids numerical explosions
  norm_factor = max(x,/nan)
  
  ;Calculate the log of the sum
  return, alog(total(exp(x-norm_factor),/nan))+norm_factor
  
End