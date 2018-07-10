      function BI1(x)
      implicit double precision (a-h,o-z)
C Abramowitz
C I1(x), -3.75<x<3.75      
      parameter (a0=0.5d0,a2=0.87890594d0,a4=.51498869d0,a6=.15084934d0,
     >a8=0.02658733d0, a10=0.00301532d0, a12=0.00032411d0)

C I1
      if (x .gt. 3.75) stop ' sbr BI1 : x .gt. 3.75'
      t = x / 3.75d0
      t2 = t * t
      BI1 = x * (a0 + (a2 + (a4 + (a6 + (a8 +  (a10 + a12 
     >                            * t2)*t2)*t2)*t2)*t2)*t2) 
      return

C I1/x
      entry x1BI1(x)
      if (x .gt. 3.75) stop ' sbr BI1 : x .gt. 3.75'
      t = x / 3.75d0
      t2 = t * t
      x1BI1 = a0 + (a2 + (a4 + (a6 + (a8 +  (a10 + a12 
     >                            * t2)*t2)*t2)*t2)*t2)*t2
      return
     
C d(I1/x) / dx, en vue de I1' = u'x/3.75 + I1/x et u'=d(I1/x) / dx = dx1I1
      entry dx1BI1(x)
      if (x .gt. 3.75) stop ' sbr BI1 : x .gt. 3.75'
      t = x / 3.75d0
      t2 = t * t
      dx1BI1=
     > (2.d0*a2 +(4.d0*a4 +(6.d0*a6 +(8.d0*a8 +(10.d0*a10 +12.d0*a12 
     >                            * t2)*t2)*t2)*t2)*t2)*t
      return     
      end
