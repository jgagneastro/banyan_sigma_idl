Function matrix_vector_product, A, v
  ;+
  ; NAME:
  ;       MATRIX_VECTOR_PRODUCT
  ;
  ; PURPOSE:
  ;       This function performs the matrix-vector multiplication of many pairs of DxD matrices or D vectors
  ;       that is optimized to perform this operation on a large number N of A matrices and v vectors at once.
  ;       This is done by looping over the dimensions D rather than the number N of matrices and vectors.
  ;
  ; INPUTS:
  ;       A - Set of N matrices of dimension DxD. Therefore A has dimensions NxDxD.
  ;       v - Set of N vectors of dimension D. Therefore v has dimensions NxD.
  ;
  ; OUTPUTS:
  ;       A set of N vectors (dimension NxD) that result from the matrix-vector multiplication
  ;       of each of the N matrices and vectors stored in A and v, respectively.
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
  if (size(v))[0] ne 2L then $
    message, ' vector must have 2 dimensions : NOBJ, NEL'
  if (size(A))[3] ne nel or (size(v))[2] ne nel then $
    message, ' Dimensions of matrix and vector do not agree !'
  nobj = (size(A))[1]

  ;w = dblarr(nobj,nel)
  w = make_array(nobj,nel,/double)
  for i=0L, nel-1L do $
    for k=0L, nel-1L do $
      w[*,i] += A[*,i,k] * v[*,k]

  return, w
End