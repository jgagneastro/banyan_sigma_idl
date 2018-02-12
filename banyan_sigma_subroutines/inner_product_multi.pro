Function inner_product_multi, u, v, metric
  ;+
  ; NAME:
  ;       INNER_PRODUCT_MULTI
  ;
  ; PURPOSE:
  ;       This function performs a vector scalar product of many pairs of D-dimension u and v vectors
  ;       induced by a single, constant metric matrix of dimension DxD. This function is optimized
  ;       to operate on a large number N of u and v vectors at once. This is done by looping over the
  ;       dimensions D rather than the number N of vectors.
  ;
  ; INPUTS:
  ;       u - Set of N vectors of dimension D. Therefore u has dimensions NxD.
  ;       v - Set of N vectors of dimension D. Therefore v has dimensions NxD.
  ;       metric - DxD matrix that induces the scalar product.
  ;
  ; OUTPUTS:
  ;       A set of N scalars (dimension N) that result from the scalar product
  ;       of each of the N u and v vectors, induced by the metric matrix.
  ;
  ; PROCEDURES USED:
  ;       NONE
  ;
  ; MODIFICATION HISTORY:
  ;       WRITTEN, Jonathan Gagne, October, 25 2017
  ;-
  
  nel = (size(u))[2]
  nobj = (size(u))[1]

  w = dblarr(nobj)
  for i=0L, nel-1L do $
    for j=0L, nel-1L do $
      w += u[*,i] * v[*,j] * metric[i,j]
  
  return, w
End