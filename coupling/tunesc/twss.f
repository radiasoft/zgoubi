c     twiss.f
c     written by Frédéric Desforges, frederic.desforges@grenoble-inp.org

      SUBROUTINE TWSS(P, 
     >                  F0,rPARAM,C)

      IMPLICIT DOUBLE PRECISION (A-H,M-Z)
      DIMENSION F0(6,6)
      DIMENSION P(4,4),C(2,2),P12(2,2),P22INV(2,2),WORK(2)
      INTEGER   IPIV,INFO
      DIMENSION IPIV(4)

      F0(1,2) =  P(2,1)/P(2,2)         ! -alpha
      F0(2,1) = F0(1,2)
      F0(1,1) = P(1,1)/P(2,2)
      F0(2,2) = (P(2,1)**2+P(2,2)**2)/(P(1,1)*P(2,2))

      F0(3,4) =  P(4,3)/P(4,4)         !  -alpha
      F0(4,3) = F0(3,4) 
      F0(3,3) = P(3,3)/P(4,4)
      F0(4,4) = (P(4,3)**2+P(4,4)**2)/(P(3,3)*P(4,4))

      rPARAM  = (SQRT(P(1,1)*P(2,2))+SQRT(P(3,3)*P(4,4)))/2

      P12(1,1) = P(1,3)
      P12(1,2) = P(1,4)
      P12(2,1) = P(2,3)
      P12(2,2) = P(2,4)
      P22INV(1,1) = P(3,3)
      P22INV(1,2) = P(3,4)
      P22INV(2,1) = P(4,3)
      P22INV(2,2) = P(4,4)
      CALL DGETRF(2,2,P22INV,2,IPIV,INFO)
      CALL DGETRI(2,P22INV,2,IPIV,WORK,2,INFO)

      C = rPARAM*MATMUL(P12,P22INV)

      RETURN                 
      END
