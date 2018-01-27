Function matrix_multiply_square_act, A, B
  ;+
  ; NAME:
  ;       MATRIX_MULTIPLY_SQUARE_ACT
  ;
  ; PURPOSE:
  ;       This function performs the matrix multiplication of many pairs of two DxD matrices that is optimized
  ;       to perform this operation on a large number of B matrices at once, and A is a single constant matrix.
  ;       This is done by looping over the dimensions D rather than the number of B matrices N.
  ;
  ; INPUTS:
  ;       A - Matrix of dimension DxD
  ;       B - Set of N matrices of dimension DxD. Therefore B has dimensions NxDxD.
  ;
  ; OUTPUTS:
  ;       A set of N matrices (dimension NxDxD) that result from the matrix multiplication.
  ;       of matrix A with each of the N matrices stored in B. .
  ;
  ; PROCEDURES USED:
  ;       NONE
  ;
  ; MODIFICATION HISTORY:
  ;       WRITTEN, Jonathan Gagne, October, 25 2017
  ;-

  if (size(B))[0] ne 3L then $
    message, ' Matrix B must have at least 3 dimensions : NOBJ, NEL, NEL'
  nel = (size(B))[2]
  if (size(A))[2] ne nel or (size(B))[3] ne nel or (size(A))[1] ne nel then $
    message, ' Dimensions of matrices do not agree !'
  nobj = (size(B))[1]

  C = dblarr(nobj,nel,nel)
  for i=0L, nel-1L do $
    for j=0L, nel-1L do $
      for k=0L, nel-1L do $
        C[*,i,j] += A[i,k] * B[*,k,j]

  return, C
End