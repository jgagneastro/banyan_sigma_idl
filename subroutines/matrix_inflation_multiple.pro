Function matrix_inflation_multiple, MATRIX, INFLATION
  ;+
  ; NAME:
  ;       MATRIX_INFLATION_MULTIPLE
  ;
  ; PURPOSE:
  ;       This function inflates the diagonal elements of a single covarianc matrix MATRIX of dimension DxD, 
  ;       using a set of N inflation vectors INFLATION of dimension D, while preserving the (off-diagonal) covariances.
  ;       This is done with the set of matrix operations INFLATION # (MATRIX # INFLATION). This function is optimized to
  ;       operate on a large number N of inflation vectors at once, by looping over the dimensions D rather than the
  ;       number N of vectors.
  ;
  ; INPUTS:
  ;       MATRIX - Matrix of dimension DxD.
  ;       INFLATION - Set of N vectors of dimension D. Therefore INFLATION has dimensions NxD.
  ;
  ; OUTPUTS:
  ;       A set of N matrices (dimension NxDxD) that result from the inflation of the diagonal
  ;       elements of MATRIX with each of the N inflation vectors INFLATION.
  ;
  ; PROCEDURES USED:
  ;       NONE
  ;
  ; MODIFICATION HISTORY:
  ;       WRITTEN, Jonathan Gagne, October, 25 2017
  ;-
  ;
  
  if (size(MATRIX))[0] ne 2L then $
    message, ' MATRIX must have 2 dimensions : NEL, NEL'
  nel = (size(MATRIX))[1]
  if (size(INFLATION))[0] ne 2L then $
    message, ' INFLATION must have 2 dimensions : NOBJ, NEL'
  nobj = (size(INFLATION))[1]
  if (size(MATRIX))[2] ne nel or (size(INFLATION))[2] ne nel then $
    message, ' Dimensions of matrices do not agree !'
  
  ;C = MATRIX # INFLATION
  C = dblarr(nobj,nel,nel)
  for i=0L, nel-1L do $
    for j=0L, nel-1L do $
      C[*,i,j] = MATRIX[i,j] * INFLATION[*,j]
  
  ;D = INFLATION # C
  D = dblarr(nobj,nel,nel)
  for i=0L, nel-1L do $
    for j=0L, nel-1L do $
      D[*,i,j] = INFLATION[*,i] * C[*,i,j]
  
  return, D
End