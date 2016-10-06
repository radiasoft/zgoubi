c     hamilt.f
c
c     written by Frédéric Desforges, frederic.desforges@grenoble-inp.org
c
c
      SUBROUTINE HAMILT(NU1,NU2,rPARAM,CMOINS,DELTA,NUX0,NUY0,DELTA2
     >,CPLUS)

      IMPLICIT DOUBLE PRECISION (A-H,M-Z)
      DOUBLE PRECISION NU1,NU2,NUX0,NUY0


      IF(rPARAM .LE. 1.D0 .AND. rPARAM .GE. SQRT(2.D0)/2.D0) THEN

            CMOINS = 2.D0*SQRT(ABS(1.D0-1.D0/rPARAM**2))/(1.D0+
     >  ABS(1.D0-1.D0/rPARAM**2))*ABS(NU1-NU2)
            DELTA  = (1.D0-ABS(1.D0-1.D0/rPARAM**2))/
     >  (1.D0+ABS(1.D0-1.D0/rPARAM**2))*(NU1-NU2)
            NUX0   = NU1 + 0.5D0*DELTA - SIGN(1.D0,DELTA)*0.5D0*
     >  SQRT(DELTA**2+CMOINS**2)
            NUY0   = NU2 - 0.5D0*DELTA + SIGN(1.D0,DELTA)*0.5D0
     >  *SQRT(DELTA**2+CMOINS**2)
            DELTA2 = (NUX0 + NUY0) - NINT(NUX0 + NUY0)
        CPLUS  = SQRT((NUY0+NUX0-INT(NUX0+NUY0))**2-(NU1+NU2-INT(NU1
     >  +NU2))**2)


      ELSE IF(rPARAM .LT. SQRT(2.)/2.) THEN

            rPARAM = rPARAM+2*(SQRT(2.D0)/2.D0-rPARAM)                        ! According to the chosen convention 
            CMOINS = 2.D0*SQRT(ABS(1.D0-1.D0/rPARAM**2))/(1.D0+               ! rPARAM cannot be lower than sqrt(2)/2
     >  ABS(1.D0-1.D0/rPARAM**2))*ABS(NU2-NU1)                                  ! For more detailed arguments look my internal
            DELTA  = (1.D0-ABS(1.D0-1.D0/rPARAM**2))/(1.D0+                   ! BNL report.
     >  ABS(1.D0-1.D0/rPARAM**2))*(NU2-NU1) 
            NUX0   = NU2 + 0.5D0*DELTA - SIGN(1.D0,DELTA)*0.5D0*
     >  SQRT(DELTA**2+CMOINS**2)
            NUY0   = NU1 - 0.5D0*DELTA + SIGN(1.D0,DELTA)*0.5D0*
     >  SQRT(DELTA**2+CMOINS**2)
            DELTA2 = (NUX0 + NUY0) - NINT(NUX0 + NUY0)
        CPLUS  = SQRT((NUY0+NUX0-INT(NUX0+NUY0))**2-(NU1+NU2-INT(NU1
     >  +NU2))**2)

      
      ELSE

            CMOINS = 0.D0
            CPLUS  = 0.D0
            NUX0   = NU1
            NUY0   = NU2
            DELTA  = 0.D0
            DELTA2 = 0.D0
c            WRITE (*,FMT='(/,''WARNING: CLOSER TO A LINEAR SUM RESONANCE
c     > THAN A LINEAR DIFFERENCE RESONANCE -> HAMILTONIAN PERTURBATION TH 
c     >EORY NOT VALID'')')

      ENDIF

      RETURN      
      END
