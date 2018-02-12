Function matrix_vector_product_vct, A, v
  ;+
  ; NAME:
  ;       MATRIX_VECTOR_PRODUCT_VCT
  ;
  ; PURPOSE:
  ;       This function performs the matrix-vector multiplication of many pairs of DxD matrices or D vectors
  ;       that is optimized to perform this operation on a large number of A matrices at once, and v is a single
  ;       constant vector. This is done by looping over the dimensions D rather than the number of A matrices N.
  ;
  ; INPUTS:
  ;       A - Set of N matrices of dimension DxD. Therefore A has dimensions NxDxD.
  ;       v - A single vector of dimension D.
  ;
  ; OUTPUTS:
  ;       A set of N vectors (dimension NxD) that result from the matrix-vector multiplication.
  ;       of each of the N matrices stored in A with vector v.
  ;
  ; PROCEDURES USED:
  ;       NONE
  ;
  ; MODIFICATION HISTORY:
  ;       WRITTEN, Jonathan Gagne, October, 25 2017
  ;-
  
  if (size(A))[0] ne 3L then $
    message, ' Matrix must have 3 dimensions : NOBJ, NEL, NEL'
  nel = (size(A))[2]
  if (size(v))[0] ne 1L then $
    message, ' vector must have 1 dimension : NEL'
  if (size(A))[3] ne nel or (size(v))[1] ne nel then $
    message, ' Dimensions of matrix and vector do not agree !'
  nobj = (size(A))[1]

  w = dblarr(nobj,nel)
  for i=0L, nel-1L do $
    for k=0L, nel-1L do $
      w[*,i] += A[*,i,k] * v[k]
  
  return, w
End