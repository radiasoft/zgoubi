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
      SUBROUTINE PROPAG(N_ET,NTRANS,NRES,NTWISS,P,C,NU1,NU2,CMOINS,
     >NUX0,NUY0)

      IMPLICIT DOUBLE PRECISION (A-H,M-Z)
      DOUBLE PRECISION G(4,4),R(4,4),P(4,4),NU1,NU2,R11(2,2),R22INV(2,2)
     >,C(2,2),WORK(2)
      INTEGER I,J,N_ET,NTRANS,NRES,NTWISS,IPIV,INFO,DEBSTR,FINSTR
      DIMENSION IPIV(4)
      CHARACTER*300 BUFFER,LABEL,KEYWOR
      LOGICAL STRCON

 1    CONTINUE
      
      READ(NTRANS,FMT='(a)',END=99,ERR=98) BUFFER
      READ(NTRANS,FMT='(a)',END=99,ERR=98) BUFFER
      READ(NTRANS,FMT='(a)',END=99,ERR=98) BUFFER
      READ(NTRANS,FMT='(a)',END=99,ERR=98) BUFFER
      READ(NTRANS,*,END=99,ERR=98) ((R(I,J),J=1,4),I=1,4)


      G = MATMUL(R,P)

c      WRITE(*,FMT='(/,6X,''MATRIX (Gij):'',/)')
c      WRITE(*,200) ((G(I,J),J=1,4),I=1,4)
c 200     FORMAT(6X,4F13.8)


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

      IF(rPARAM .LE. 1. .AND. rPARAM .GE. SQRT(2.)/2.) THEN
            CMOINS = 2*SQRT(ABS(1-1/rPARAM**2))/(1+ABS(1-1/rPARAM**2))*A
     >BS(NU1-NU2)
      ELSE IF(rPARAM .LT. SQRT(2.)/2.) THEN
            rPARAM = rPARAM+2*(SQRT(2.)/2.-rPARAM)                        ! According to the chosen convention 
            CMOINS = 2*SQRT(ABS(1-1/rPARAM**2))/(1+ABS(1-1/rPARAM**2))*A  ! rPARAM cannot be lower than sqrt(2)/2
     >BS(NU2-NU1)                                                         ! For more detailed arguments look my internal
      ELSE
            CMOINS = 0.
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

c      WRITE(*,FMT='(/,6X,''rPARAM, BETA1, BETA2:'',3F15.8)') rPARAM,BETA
c     >1,BETA2


c     Reading of zgoubi.res for the position corresponding to the transported Twiss' parameters 

      CALL POSITI(NRES,
     >                 ARCLEN,LABEL,KEYWOR)


c     Extraction of the generalized Twiss' parameters at different positions

      CALL EXTRAC(N_ET,ARCLEN,R,rPARAM,C,NU1,NU2,ALPHA1,ALPHA2,BETA1,BET
     >A2,GAMMA1,GAMMA2,CMOINS,CPLUS,DELTA,DELTA2,NUX0,NUY0,P)
      
c      WRITE(NTWISS,FMT='(8(F13.8,X))') ARCLEN/100,BETA1,BETA2,ALPHA1,ALP
c     >HA2,GAMMA1,GAMMA2,rPARAM
      WRITE(NTWISS,FMT='((a),6X,(a),8(F13.8,X))') LABEL(1:9),KEYWOR(1:20
     >),ARCLEN/100,BETA1,BETA2,ALPHA1,ALPHA2,GAMMA1,GAMMA2,rPARAM

      GOTO 1


 98   CONTINUE
      
      WRITE(*,FMT='(/,/,''ERROR IN TUNESC.F DURING THE READING OF TRANSF
     >ERM.DAT'',/,/)')

 99   CONTINUE

       
      END
