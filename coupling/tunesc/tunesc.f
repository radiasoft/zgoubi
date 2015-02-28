C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  François Méot
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.
C
C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 51 Franklin Street, Fifth Floor,
C  Boston, MA  02110-1301  USA
C
C  François Meot <fmeot@bnl.gov>
C  Brookhaven National Laboratory 
C  C-AD, Bldg 911
C  Upton, NY, 11973, USA
C  -------
c     tunesc.f
c     written by Frédéric Desforges, frederic.desforges@grenoble-inp.org
c
C     Computation of the eigentunes, the generalized Twiss' parameters
c     and the perturbation hamiltonian theory parameters in case of 
c     weak linear difference coupling resonances
c
C     Assumptions:
c                 - symplecticity
c                 - weak and linear coupling 
C                 - integer and half-integer tunes forbidden
C                 - vertical tune and horizontal not perfectly equal
C
c
      SUBROUTINE TUNESC(R6,
     >                     F0,NU1,NU2,CMUY,CMUZ,IERY,IERZ,RPARAM,C)
      IMPLICIT DOUBLE PRECISION (A-H,M-Z)
      DIMENSION F0(6,6), R6(6,6)
      DIMENSION R(4,4),REIG(4,4),Q(4,4),U(4,4),WR(4),WI(4),P(4,4),PINV(4
     >,4),A(4,4),GV1(16),FV1(16),WORK(4),G(4,4),C(2,2),NUTEST(4)
      DOUBLE COMPLEX EV1,EV2,EV3,EV4,V(4,4),PHASE1,PHASE2,F1001,F1010
      INTEGER   IPIV,IERR,INFO,IOS1,IOS2,IOS3,IOS4,IOS7,IOS8
      DIMENSION IPIV(4)
      CHARACTER(300) BUFFER
      LOGICAL STRCON,N_PROP
      CHARACTER(300) LABEL,KEYWOR

      PARAMETER (PI = 3.141592653589793)
      logical ok, idluni

      IERY=0
      IERZ=0

      CALL RAZ(F0, 6*6)
      DO J = 1, 4
        DO I = 1, 4
          R(I,J) = R6(I,J)
        ENDDO
      ENDDO

C     Computation of eigen values and eigen vectors
      
      REIG = R
      CALL RG(4,4,REIG,WR,WI,1,Q,GV1,FV1,IERR)         ! The output matrix Q is composed of as follows:
                                                       ! 1st column: real part of eigenvector 1
      EV1 = CMPLX(WR(1),WI(1))                       ! 2nd column: imaginary part of eigenvector 1
      EV2 = CMPLX(WR(2),WI(2))                       ! 3rd column: real part of eigenvector 2
      EV3 = CMPLX(WR(3),WI(3))                       ! 4th column: imaginary part of eigenvector 2
      EV4 = CMPLX(WR(4),WI(4))

C     Normalization of the eigen vectors

      NORM1  = SQRT(Q(1,1)**2+Q(1,2)**2+Q(2,1)**2+Q(2,2)**2+Q(3,1)**2+Q(
     >3,2)**2+Q(4,1)**2+Q(4,2)**2)
      NORM2  = SQRT(Q(1,3)**2+Q(1,4)**2+Q(2,3)**2+Q(2,4)**2+Q(3,3)**2+Q(
     >3,4)**2+Q(4,3)**2+Q(4,4)**2)

      V(1,1) = CMPLX(Q(1,1),Q(1,2))/NORM1
      V(2,1) = CMPLX(Q(2,1),Q(2,2))/NORM1
      V(3,1) = CMPLX(Q(3,1),Q(3,2))/NORM1
      V(4,1) = CMPLX(Q(4,1),Q(4,2))/NORM1
      V(1,2) = CMPLX(Q(1,1),-Q(1,2))/NORM1
      V(2,2) = CMPLX(Q(2,1),-Q(2,2))/NORM1
      V(3,2) = CMPLX(Q(3,1),-Q(3,2))/NORM1
      V(4,2) = CMPLX(Q(4,1),-Q(4,2))/NORM1
      V(1,3) = CMPLX(Q(1,3),Q(1,4))/NORM2
      V(2,3) = CMPLX(Q(2,3),Q(2,4))/NORM2
      V(3,3) = CMPLX(Q(3,3),Q(3,4))/NORM2
      V(4,3) = CMPLX(Q(4,3),Q(4,4))/NORM2
      V(1,4) = CMPLX(Q(1,3),-Q(1,4))/NORM2
      V(2,4) = CMPLX(Q(2,3),-Q(2,4))/NORM2
      V(3,4) = CMPLX(Q(3,3),-Q(3,4))/NORM2
      V(4,4) = CMPLX(Q(4,3),-Q(4,4))/NORM2


C     Phase uncertainity

      THETA1 = ATAN2(AIMAG(V(1,1)),REAL(V(1,1)))
      THETA2 = ATAN2(AIMAG(V(3,3)),REAL(V(3,3)))
      PHASE1 = CMPLX(COS(-PI/2.D0-THETA1),SIN(-PI/2.D0-THETA1))
      PHASE2 = CMPLX(COS(-PI/2.D0-THETA2),SIN(-PI/2.D0-THETA2))

      V(1,1) = CMPLX(Q(1,1),Q(1,2))/NORM1*PHASE1
      V(2,1) = CMPLX(Q(2,1),Q(2,2))/NORM1*PHASE1
      V(3,1) = CMPLX(Q(3,1),Q(3,2))/NORM1*PHASE1
      V(4,1) = CMPLX(Q(4,1),Q(4,2))/NORM1*PHASE1
      V(1,2) = CMPLX(Q(1,1),-Q(1,2))/NORM1*CONJG(PHASE1)
      V(2,2) = CMPLX(Q(2,1),-Q(2,2))/NORM1*CONJG(PHASE1)
      V(3,2) = CMPLX(Q(3,1),-Q(3,2))/NORM1*CONJG(PHASE1)
      V(4,2) = CMPLX(Q(4,1),-Q(4,2))/NORM1*CONJG(PHASE1)
      V(1,3) = CMPLX(Q(1,3),Q(1,4))/NORM2*PHASE2
      V(2,3) = CMPLX(Q(2,3),Q(2,4))/NORM2*PHASE2
      V(3,3) = CMPLX(Q(3,3),Q(3,4))/NORM2*PHASE2
      V(4,3) = CMPLX(Q(4,3),Q(4,4))/NORM2*PHASE2
      V(1,4) = CMPLX(Q(1,3),-Q(1,4))/NORM2*CONJG(PHASE2)
      V(2,4) = CMPLX(Q(2,3),-Q(2,4))/NORM2*CONJG(PHASE2)
      V(3,4) = CMPLX(Q(3,3),-Q(3,4))/NORM2*CONJG(PHASE2)
      V(4,4) = CMPLX(Q(4,3),-Q(4,4))/NORM2*CONJG(PHASE2)

 
C      Computation of the fractional part of the horizontal tunes without considering coupling

      CMUY = .5D0 * (R(1,1)+R(2,2))
      SMUY = SIGN(SQRT(-R(1,2)*R(2,1)-.25D0*(R(1,1)-R(2,2))**2),R(1,2))
      NUY  = SIGN(ATAN2(SMUY,CMUY) /(2.D0 * PI) ,R(1,2))
      IF(NUY .LT. 0.D0) NUY = 1 + NUY

C     Computation of the fractional part of the vertical tunes without considering coupling

      CMUZ = .5D0 * (R(3,3)+R(4,4))
      SMUZ = SIGN(SQRT(-R(3,4)*R(4,3)-.25D0*(R(3,3)-R(4,4))**2),R(3,4))
      NUZ  = SIGN(ATAN2(SMUZ,CMUZ) /(2.D0 * PI) ,R(3,4))
      IF(NUZ .LT. 0.D0) NUZ = 1 + NUZ

C     Tunes in the decoupled referential

      NUTEST(1) = ATAN2(WI(1),WR(1))/(2.D0*PI)          ! It is important to notice that with this 
            IF (NUTEST(1) .LT. 0.0D0) THEN              ! numerical process one imposes NU1 to be 
                NUTEST(1) = - NUTEST(1)                 ! perturbed NUX and NU2 the perturbed NY
                NUTEST(2) = 1 - NUTEST(1)
            ELSE
                NUTEST(2) = 1 - NUTEST(1)
            ENDIF
      NUTEST(3) = ATAN2(WI(3),WR(3))/(2.D0*PI)
            IF (NUTEST(2) .LT. 0.0D0) THEN
                NUTEST(3) = - NUTEST(3)
                NUTEST(4) = 1 - NUTEST(3)
            ELSE
                NUTEST(4) = 1 - NUTEST(3)
            ENDIF

      NU1 = NUTEST(1)
      D1  = ABS(NUTEST(1)-NUY)
      IF(ABS(NUTEST(1)-NUZ) .LT. D1) D1=ABS(NUTEST(1)-NUZ) 
      D2  = ABS(NUTEST(2)-NUY)
      IF(ABS(NUTEST(2)-NUZ) .LT. D2) D2=ABS(NUTEST(2)-NUZ)
      IF (D2 .LT. D1) NU1 = NUTEST(2)

      NU2 = NUTEST(3)
      D3  = ABS(NUTEST(3)-NUY)
      IF(ABS(NUTEST(3)-NUZ) .LT. D3) D3=ABS(NUTEST(3)-NUZ) 
      D4  = ABS(NUTEST(4)-NUY)
      IF(ABS(NUTEST(4)-NUZ) .LT. D4) D4=ABS(NUTEST(4)-NUZ)
      IF (D4 .LT. D3) NU2 = NUTEST(4)


C     Transfert matrix in action-angle referential

      COSMU1 = (WR(1)+WR(2))/2
      SINMU1 = (WI(1)-WI(2))/2
      IF(NU1 .EQ. NUTEST(2)) SINMU1 = -SINMU1

      COSMU2 = (WR(3)+WR(4))/2
      SINMU2 = (WI(3)-WI(4))/2
      IF(NU2 .EQ. NUTEST(4)) SINMU2 = -SINMU2

      U(1,1) = COSMU1
      U(1,2) = SINMU1
      U(1,3) = 0.D0
      U(1,4) = 0.D0
      U(2,1) = -SINMU1
      U(2,2) = COSMU1
      U(2,3) = 0.D0
      U(2,4) = 0.D0
      U(3,1) = 0.D0
      U(3,2) = 0.D0
      U(3,3) = COSMU2
      U(3,4) = SINMU2
      U(4,1) = 0.D0
      U(4,2) = 0.D0
      U(4,3) = -SINMU2
      U(4,4) = COSMU2

c      WRITE(*, FMT='(/,/,6X,''TRANSFERT MATRIX (Uij) IN THE ANGLE-ACTIO
c     >N FRAME'',/)')
c      WRITE(*,200) (( U(I,J) , J=1,4) , I=1,4)
c 200     FORMAT(6X,4F13.8)


C     Transformation matrix from coupled referential to action-angle one

      P(1,1) = DBLE(CMPLX(0.D0,1.D0/SQRT(2.D0))*(V(1,1)-V(1,2)))
      P(2,1) = DBLE(CMPLX(0.D0,1.D0/SQRT(2.D0))*(V(2,1)-V(2,2)))
      P(3,1) = DBLE(CMPLX(0.D0,1.D0/SQRT(2.D0))*(V(3,1)-V(3,2)))
      P(4,1) = DBLE(CMPLX(0.D0,1.D0/SQRT(2.D0))*(V(4,1)-V(4,2)))
      P(1,2) = DBLE(CMPLX(1.D0/SQRT(2.D0),0.D0)*(V(1,1)+V(1,2)))
      P(2,2) = DBLE(CMPLX(1.D0/SQRT(2.D0),0.D0)*(V(2,1)+V(2,2)))
      P(3,2) = DBLE(CMPLX(1.D0/SQRT(2.D0),0.D0)*(V(3,1)+V(3,2)))
      P(4,2) = DBLE(CMPLX(1.D0/SQRT(2.D0),0.D0)*(V(4,1)+V(4,2)))
      P(1,3) = DBLE(CMPLX(0.D0,1.D0/SQRT(2.D0))*(V(1,3)-V(1,4)))
      P(2,3) = DBLE(CMPLX(0.D0,1.D0/SQRT(2.D0))*(V(2,3)-V(2,4)))
      P(3,3) = DBLE(CMPLX(0.D0,1.D0/SQRT(2.D0))*(V(3,3)-V(3,4)))
      P(4,3) = DBLE(CMPLX(0.D0,1.D0/SQRT(2.D0))*(V(4,3)-V(4,4)))
      P(1,4) = DBLE(CMPLX(1.D0/SQRT(2.D0),0.D0)*(V(1,3)+V(1,4)))
      P(2,4) = DBLE(CMPLX(1.D0/SQRT(2.D0),0.D0)*(V(2,3)+V(2,4)))
      P(3,4) = DBLE(CMPLX(1.D0/SQRT(2.D0),0.D0)*(V(3,3)+V(3,4)))
      P(4,4) = DBLE(CMPLX(1.D0/SQRT(2.D0),0.D0)*(V(4,3)+V(4,4)))

      IF(NU1 .EQ. NUTEST(2)) P(:,1) = -P(:,1)
      IF(NU2 .EQ. NUTEST(4)) P(:,3) = -P(:,3)

      CALL GETDET(P,4,DETP)

      P=P/DETP**(0.25)

c      WRITE(*,FMT='(/,/,6X,''TRANSFORMATION MATRIX (Pij)'',/)')
c      WRITE(*,300) (( P(I,J) , J=1,4) , I=1,4)
c 300     FORMAT(6X,4F13.8)
    
 
C     Inverse of the transformation matrix

      PINV = P
      CALL DGETRF(4,4,PINV,4,IPIV,INFO)
      CALL DGETRI(4,PINV,4,IPIV,WORK,4,INFO)


C     Edwards-Teng's parameters: alpha, beta, gamma, r and the matrix C

      CALL TWSS(P, 
     >            F0,rPARAM,C)

C     Hamiltonian pertubation parameters and  Computing of the unperturbed tunes from the coupling parameters
      
      CALL HAMILT(NU1,NU2,rPARAM,CMOINS,DELTA,NUX0,NUY0,DELTA2,CPLUS)

       call matic2(rPARAM,C,NU1,NU2,ALPHA1,ALPHA2,BETA1,BETA2,
     > GAMMA1,GAMMA2,CMOINS,CPLUS,DELTA,DELTA2,NUX0,NUY0,P)      

C     Propagation of the generalized Twiss' parameters and coupling parameters

      N_PROP = .FALSE.                                            
C      N_PROP = .TRUE.                                                    ! Flag turned on by default,

      IF(N_PROP .EQV. .TRUE.) THEN

        OPEN(3, FILE='ETbeta.res',STATUS='UNKNOWN',IOSTAT=IOS3)
        OPEN(4, FILE='transfertM.dat',STATUS='UNKNOWN',IOSTAT=IOS4)
        OPEN(7, FILE='zgoubi.res',STATUS='UNKNOWN',IOSTAT=IOS7)
        OPEN(8, FILE='twiss.out',STATUS='UNKNOWN',IOSTAT=IOS8)
      
        WRITE(8,FMT='(2X,''LABEL'',9X,''KEYWORD'',20X,''S'',10X,
     >  ''BETA1'',9X,''BETA2'',9X,''ALPHA1'',8X,''ALPHA2'',8X,
     >  ''GAMMA1'',8X,''GAMMA2'',8X,''rPARAM'',/)')
      
        CALL PROPAG(3,4,nres,8,P,C,NU1,NU2,CMOINS,NUX0,NUY0)    ! the generalized Twiss' parameters 
                                                                         ! and coupling strenght along 
        CLOSE(8,IOSTAT=IOS8)                                               ! the ring
        CLOSE(7,IOSTAT=IOS7)      
        CLOSE(4,IOSTAT=IOS4)
        CLOSE(3,IOSTAT=IOS3)
      ENDIF

      RETURN
      END
