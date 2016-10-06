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

      SUBROUTINE TUNESC(RT6,
     >                     F0,NU1,NU2,CMUY,CMUZ,IERY,IERZ)
      IMPLICIT DOUBLE PRECISION (A-H,M-Z)
      DIMENSION RT6(6,6)
      DIMENSION F0(6,6)
      DIMENSION Q(4,4),U(4,4),WR(4),WI(4)
     >,PINV(4,4),GV1(16),FV1(16),WORK(4),NUTEST(4)
C FM Apr. 2015
      INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15,50)
      COMPLEX(DP) EV1,EV2,EV3,EV4,V(4,4),PHASE1,PHASE2
      INTEGER   IPIV,IERR,INFO
      DIMENSION IPIV(4)

      DIMENSION RT(4,4), REIG(4,4)
      DIMENSION P(4,4),C(2,2)
      DIMENSION PO(4,4),CO(2,2)
      SAVE P, C

      PARAMETER (PI = 3.141592653589793D0)

      IERY=0
      IERZ=0

      CALL RAZ(F0, 6*6)

C RT6 is the 6*6 1-turn map
      DO J = 1, 4
        DO I = 1, 4
          RT(I,J) = RT6(I,J)
        ENDDO
      ENDDO
C     Computation of eigen values and eigen vectors

      REIG = RT
      CALL RG(4,4,REIG,WR,WI,1,Q,GV1,FV1,IERR)         ! The output matrix Q is composed of as follows:
                                                     ! 1st column: real part of eigenvector 1
      EV1 = CMPLX(WR(1),WI(1))                       ! 2nd column: imaginary part of eigenvector 1
      EV2 = CMPLX(WR(2),WI(2))                       ! 3rd column: real part of eigenvector 2
      EV3 = CMPLX(WR(3),WI(3))                       ! 4th column: imaginary part of eigenvector 2
      EV4 = CMPLX(WR(4),WI(4))

C     Normalization of the eigen vectors

      XNRM1  = SQRT(Q(1,1)**2+Q(1,2)**2+Q(2,1)**2+Q(2,2)**2+Q(3,1)**2+Q(
     >3,2)**2+Q(4,1)**2+Q(4,2)**2)
      XNRM2  = SQRT(Q(1,3)**2+Q(1,4)**2+Q(2,3)**2+Q(2,4)**2+Q(3,3)**2+Q(
     >3,4)**2+Q(4,3)**2+Q(4,4)**2)

C             write(*,*) ' Q11, Q12 ',Q(1,1),Q(1,2),XNRM1
C             write(*,*) ' CMPLX(Q(1,1),Q(1,2)) ',CMPLX(Q(1,1),Q(1,2))
C             write(*,*) '  ',XNRM1,CMPLX(Q(1,1),Q(1,2))/XNRM1

      V(1,1) = CMPLX(Q(1,1),Q(1,2))/XNRM1
      V(2,1) = CMPLX(Q(2,1),Q(2,2))/XNRM1
      V(3,1) = CMPLX(Q(3,1),Q(3,2))/XNRM1
      V(4,1) = CMPLX(Q(4,1),Q(4,2))/XNRM1
      V(1,2) = CMPLX(Q(1,1),-Q(1,2))/XNRM1
      V(2,2) = CMPLX(Q(2,1),-Q(2,2))/XNRM1
      V(3,2) = CMPLX(Q(3,1),-Q(3,2))/XNRM1
      V(4,2) = CMPLX(Q(4,1),-Q(4,2))/XNRM1
      V(1,3) = CMPLX(Q(1,3),Q(1,4))/XNRM2
      V(2,3) = CMPLX(Q(2,3),Q(2,4))/XNRM2
      V(3,3) = CMPLX(Q(3,3),Q(3,4))/XNRM2
      V(4,3) = CMPLX(Q(4,3),Q(4,4))/XNRM2
      V(1,4) = CMPLX(Q(1,3),-Q(1,4))/XNRM2
      V(2,4) = CMPLX(Q(2,3),-Q(2,4))/XNRM2
      V(3,4) = CMPLX(Q(3,3),-Q(3,4))/XNRM2
      V(4,4) = CMPLX(Q(4,3),-Q(4,4))/XNRM2

C     Phase uncertainity

      THETA1 = ATAN2(AIMAG(V(1,1)),REAL(V(1,1)))
      THETA2 = ATAN2(AIMAG(V(3,3)),REAL(V(3,3)))
      PHASE1 = CMPLX(COS(-PI/2.D0-THETA1),SIN(-PI/2.D0-THETA1))
      PHASE2 = CMPLX(COS(-PI/2.D0-THETA2),SIN(-PI/2.D0-THETA2))

C             write(*,*) ' v11 ',V(1,1)
C             write(*,*) ' theta1, 2 ',theta1,theta2
C             write(*,*) ' phase1, 2 ',phase1,phase2

      V(1,1) = CMPLX(Q(1,1),Q(1,2))/XNRM1*PHASE1
      V(2,1) = CMPLX(Q(2,1),Q(2,2))/XNRM1*PHASE1
      V(3,1) = CMPLX(Q(3,1),Q(3,2))/XNRM1*PHASE1
      V(4,1) = CMPLX(Q(4,1),Q(4,2))/XNRM1*PHASE1
      V(1,2) = CMPLX(Q(1,1),-Q(1,2))/XNRM1*CONJG(PHASE1)
      V(2,2) = CMPLX(Q(2,1),-Q(2,2))/XNRM1*CONJG(PHASE1)
      V(3,2) = CMPLX(Q(3,1),-Q(3,2))/XNRM1*CONJG(PHASE1)
      V(4,2) = CMPLX(Q(4,1),-Q(4,2))/XNRM1*CONJG(PHASE1)
      V(1,3) = CMPLX(Q(1,3),Q(1,4))/XNRM2*PHASE2
      V(2,3) = CMPLX(Q(2,3),Q(2,4))/XNRM2*PHASE2
      V(3,3) = CMPLX(Q(3,3),Q(3,4))/XNRM2*PHASE2
      V(4,3) = CMPLX(Q(4,3),Q(4,4))/XNRM2*PHASE2
      V(1,4) = CMPLX(Q(1,3),-Q(1,4))/XNRM2*CONJG(PHASE2)
      V(2,4) = CMPLX(Q(2,3),-Q(2,4))/XNRM2*CONJG(PHASE2)
      V(3,4) = CMPLX(Q(3,3),-Q(3,4))/XNRM2*CONJG(PHASE2)
      V(4,4) = CMPLX(Q(4,3),-Q(4,4))/XNRM2*CONJG(PHASE2)

C              write(*,*) ' tunesc ',XNRM1,PHASE1,XNRM2,PHASE2
C              write(*,*) v

C      Computation of the fractional part of the horizontal tunes without considering coupling

      CMUY = .5D0 * (RT(1,1)+RT(2,2))
      SMUY = SIGN(SQRT(-RT(1,2)*RT(2,1)
     >                    -.25D0*(RT(1,1)-RT(2,2))**2),RT(1,2))
      NUY  = SIGN(ATAN2(SMUY,CMUY) /(2.D0 * PI) ,RT(1,2))
      IF(NUY .LT. 0.D0) NUY = 1.D0 + NUY


C     Computation of the fractional part of the vertical tunes without considering coupling

      CMUZ = .5D0 * (RT(3,3)+RT(4,4))
      SMUZ = SIGN(SQRT(-RT(3,4)*RT(4,3)
     >                    -.25D0*(RT(3,3)-RT(4,4))**2),RT(3,4))
      NUZ  = SIGN(ATAN2(SMUZ,CMUZ) /(2.D0 * PI) ,RT(3,4))
      IF(NUZ .LT. 0.D0) NUZ = 1.D0 + NUZ


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

C              write(*,*) ' tunesc , p '
C              write(*,*) p
C                    read(*,*)

      CALL GETDET(P,4,DETP)

      P=P/DETP**(0.25D0)




C     Inverse of the transformation matrix

      PINV = P
      CALL DGETRF(4,4,PINV,4,IPIV,INFO)
      CALL DGETRI(4,PINV,4,IPIV,WORK,4,INFO)


C     Edwards-Teng's parameters: alpha, beta, gamma, r and the matrix C         
C      CALL TWSS(ALPHA1,BETA1,GAMMA1,ALPHA2,BETA2,GAMMA2,rPARAM,C,P)

      CALL TWSSFC(P,RT6, 
     >                 F0,rPARAM,C)
      
C     Hamiltonian pertubation parameters and  Computing of the unperturbed tunes from the coupling parameters
      
      CALL HAMILT(NU1,NU2,rPARAM,CMOINS,DELTA,NUX0,NUY0,DELTA2,CPLUS)

C For possible write out to zgoubi.res
      CALL MATIC2(RPARAM,C,NU1,NU2,-F0(1,2),-F0(3,4),F0(1,1),F0(3,3),
     >F0(2,2),F0(4,4),CMOINS,CPLUS,DELTA,DELTA2,NUX0,NUY0,P)      

      RETURN

      ENTRY TUNES1(
     >             PO,CO)
      DO J = 1, 4
        DO I = 1, 4
          PO(I,J) = P(I,J)
        ENDDO
      ENDDO
      DO J = 1, 2
        DO I = 1, 2
          CO(I,J) = C(I,J)
        ENDDO
      ENDDO
      RETURN
      END
