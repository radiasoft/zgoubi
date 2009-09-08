      function BI1p(x)
      implicit double precision (a-h,o-z)
C Abramowitz
C I1(x), -3.75<x<3.75      
C dI1/dx
      BI1p = dx1BI1(x)*x/3.75d0 + x1BI1(x)
      return
      end
