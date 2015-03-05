c     propag.f
c
c     written by Frédéric Desforges, frederic.desforges@grenoble-inp.org
c
c     The propagation of the generalized twiss parameters uses a different 
c     mathematical process from the one which was presented in my report.
c     it uses the one given by Yun Luo in the paper entitled "Linear coupling
c     parametrization in the action-angle frame" of the PHYSICAL REVIEW
c     SPECIAL TOPICS - ACCELERATORS AND BEAMS, published 7 December 2004
c
c
      SUBROUTINE PROPAG(r66,P,C,NU1,NU2,CMOINS,rparam,
     >NUX0,NUY0)

      IMPLICIT DOUBLE PRECISION (A-H,M-Z)
      dimension r66(6,6)
      DOUBLE PRECISION G(4,4),R(4,4),P(4,4),NU1,NU2,R11(2,2),R22INV(2,2)
     >,C(2,2),WORK(2)
      INTEGER I,J,N_ET,NTRANS,NRES,NTWISS,IPIV,INFO,DEBSTR,FINSTR
      DIMENSION IPIV(4)
      CHARACTER*300 BUFFER,LABEL,KEYWOR
      LOGICAL STRCON

      do j = 1, 4
        do i = 1, 4
          R(I,J) = r66(i,j)    
        enddo
      enddo

      G = MATMUL(R,P)

C     Propagation of the coupling parameters and the coupling strenght

      R11(1,1) = P(1,1)
      R11(1,2) = P(1,2)
      R11(2,1) = P(2,1)
      R11(2,2) = P(2,2)
      R22INV(1,1) = P(3,3)
      R22INV(1,2) = P(3,4)
      R22INV(2,1) = P(4,3)
      R22INV(2,2) = P(4,4)
      CALL DGETRF(2,2,R22INV,2,IPIV,INFO)
      CALL DGETRI(2,R22INV,2,IPIV,WORK,2,INFO)

      C = MATMUL(MATMUL(R11,C),R22INV)

      rPARAM = 0.5*(SQRT(G(1,1)*G(2,2)-G(1,2)*G(2,1))+SQRT(G(3,3)*G(4,4)
     >-G(3,4)*G(4,3)))

      IF(rPARAM .LE. 1.D0 .AND. rPARAM .GE. SQRT(2.D0)/2.D0) THEN
        CMOINS = 2D0*SQRT(ABS(1D0-1D0/rPARAM**2))/
     >  (1D0+ABS(1D0-1D0/rPARAM**2))*ABS(NU1-NU2)
      ELSE IF(rPARAM .LT. SQRT(2.D0)/2.D0) THEN
        rPARAM = rPARAM+2*(SQRT(2.D0)/2.D0-rPARAM)                        ! According to the chosen convention 
        CMOINS = 2.D0*SQRT(ABS(1.D0-1.D0/rPARAM**2))/
     >  (1.D0+ABS(1.D0-1.D0/rPARAM**2))*ABS(NU2-NU1) ! rPARAM cannot be lower than sqrt(2)/2     
                                                     ! For more detailed arguments look my internal
      ELSE
        CMOINS = 0.D0
      ENDIF


C     Propagation of the generalized Twiss' parameters and coupling parameters
      
      BETA1  = (G(1,1)**2+G(1,2)**2) / (G(1,1)*G(2,2)-G(1,2)*G(2,1))
      BETA2  = (G(3,3)**2+G(3,4)**2) / (G(3,3)*G(4,4)-G(3,4)*G(4,3))
      ALPHA1 = -(G(1,1)*G(2,1)+G(1,2)*G(2,2)) / (G(1,1)*G(2,2)-G(1,2)*G(
     >2,1))
      ALPHA2 = -(G(3,3)*G(4,3)+G(3,4)*G(4,4)) / (G(3,3)*G(4,4)-G(3,4)*G(
     >4,3))
      GAMMA1 = (1.+ALPHA1**2)/BETA1
      GAMMA2 = (1.+ALPHA2**2)/BETA2

      RETURN       
      END
