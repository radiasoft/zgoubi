c     twiss.f
c     written by Frédéric Desforges, frederic.desforges@grenoble-inp.org

      SUBROUTINE TWSS(ALPHA1,BETA1,GAMMA1,ALPHA2,BETA2,GAMMA2,rPARAM,C,
     >P)

      IMPLICIT DOUBLE PRECISION (A-H,M-Z)
      DIMENSION P(4,4),C(2,2),P12(2,2),P22INV(2,2),WORK(2)
      INTEGER   IPIV,INFO
      DIMENSION IPIV(4)

      ALPHA1  = -P(2,1)/P(2,2)
      BETA1   = P(1,1)/P(2,2)
      GAMMA1  = (P(2,1)**2+P(2,2)**2)/(P(1,1)*P(2,2))

      ALPHA2  = -P(4,3)/P(4,4)
      BETA2   = P(3,3)/P(4,4)
      GAMMA2  = (P(4,3)**2+P(4,4)**2)/(P(3,3)*P(4,4))

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
            
      
      END
