Function matrix_vector_error_product, A, v_error
  ;+
  ; NAME:
  ;       MATRIX_VECTOR_ERROR_PRODUCT
  ;
  ; PURPOSE:
  ;       This function performs the measurement error of the matrix-vector multiplication of many pairs of DxD
  ;       matrices or D vectors that is optimized to perform this operation on a large number N of A matrices and
  ;       v vectors at once. This is done by looping over the dimensions D rather than the number N of matrices and vectors.
  ;
  ; INPUTS:
  ;       A - Set of N matrices of dimension DxD. Therefore A has dimensions NxDxD.
  ;       v_error - Measurement errors of the set of N vectors of dimension D. Therefore v_error has dimensions NxD.
  ;
  ; OUTPUTS:
  ;       A set of N vectors (dimension NxD) that result from the measurement errors of the matrix-vector multiplication
  ;       of each of the N matrices and vectors A and v, respectively.
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
  if (size(v_error))[0] ne 2L then $
    message, ' vector must have 2 dimensions : NOBJ, NEL'
  if (size(A))[3] ne nel or (size(v_error))[2] ne nel then $
    message, ' Dimensions of matrix and vector do not agree !'
  nobj = (size(A))[1]

  ;w = dblarr(nobj,nel)
  w2 = make_array(nobj,nel,/double)
  for i=0L, nel-1L do $
    for k=0L, nel-1L do $
    w2[*,i] += (A[*,i,k] * v_error[*,k])^2d0

  return, sqrt(w2)
End