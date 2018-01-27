Function parabolic_cylinder_function, n, x, FPRIME=fprime
  ;+
  ; NAME:
  ;       PARABOLIC_CYLINDER_FUNCTION
  ;
  ; PURPOSE:
  ;       Returns the parabolic cylinder function D_n(x) with a negative integer order n
  ;       between -1 and -7. 
  ;
  ; INPUTS:
  ;       n - Order.
  ;       x - Abscissa.
  ;       
  ; OPTIONAL INPUT KEYWORD:
  ;       /FPRIME - Returns a parabolic cylinder function modified for numerical stability, by multiplying
  ;                 a factor exp(x^2/4) to its output. 
  ;
  ; PROCEDURES USED:
  ;       NONE
  ;
  ; MODIFICATION HISTORY:
  ;       WRITTEN, Jonathan Gagne, October, 25 2017
  ;-
  
  factor = keyword_set(fprime) ? 1d0 : exp(x^2/4d0)
  
  sq2 = sqrt(2d0)
  spi2 = sqrt(!dpi)/sq2
  xsq2 = x / sq2
  verfc = erfc(xsq2)
  eps = exp(-x^2/2d0)
  case n of
    -1: out = factor * spi2 * verfc
    -2: out = -1d0 * factor * (spi2*x*verfc - eps)
    -3: out = factor / 2d0 * (spi2*(x^2+1)*verfc - eps*x)
    -4: out = -1d0 * factor / 6d0 * (spi2*(x^3+3*x)*verfc - eps*(x^2+2))
    -5: out = factor / 24d0 * (spi2*(x^4+6d0*x^2+3d0)*verfc - eps*(x^3+5d0*x))
    -6: out = -1d0 * factor / 120d0 * (spi2*(x^5+10*x^3+15*x)*verfc - eps*(x^4+9*x^2+8))
    -7: out = factor / 720d0 * (spi2*(x^6+15*x^4+45*x^2+15)*verfc - eps*(x^5+14*x^3+33))
    else : message, ' This order was not built into this function !'
  endcase
  
  return, out
End