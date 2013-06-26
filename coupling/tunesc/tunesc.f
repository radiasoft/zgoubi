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
      IMPLICIT DOUBLE PRECISION (A-H,M-Z)
      DIMENSION R(4,4),REIG(4,4),Q(4,4),U(4,4),WR(4),WI(4),P(4,4),PINV(4
     >,4),A(4,4),GV1(16),FV1(16),WORK(4),G(4,4),C(2,2),NUTEST(4)
      DOUBLE COMPLEX EV1,EV2,EV3,EV4,V(4,4),PHASE1,PHASE2,F1001,F1010
      INTEGER   IPIV,IERR,INFO,IOS1,IOS2,IOS3,IOS4,IOS7,IOS8
      DIMENSION IPIV(4)
      CHARACTER*300 BUFFER
      LOGICAL STRCON,N_PROP
      CHARACTER*300 LABEL,KEYWOR

      PARAMETER (PI = 3.141592653589793)


c      CALL SYSTEM('rm transfertM.dat')
c      CALL SYSTEM('touch transfertM.dat')
      

c     Opening of the file containing the transfert matrix

      OPEN(1, FILE='transfertM.dat',STATUS='UNKNOWN',IOSTAT=IOS1)


c     Zgoubi calculations

c      CALL SYSTEM('/home/frederic/zgoubiCoupling/zgoubi/zgoubi')


 1    CONTINUE


C     Transfert matrix in coupled referential

      READ(1,FMT='(a)',END=99,ERR=98) BUFFER
      READ(1,FMT='(a)',END=99,ERR=98) BUFFER
      READ(1,FMT='(a)',END=99,ERR=98) BUFFER
      READ(1,FMT='(a)',END=99,ERR=98) BUFFER
      READ(1,*)((R(I,J),J=1,4),I=1,4)

      GOTO 1

 98   CONTINUE
      
      WRITE(*,FMT='(/,/,''ERROR IN READING OF ZGOUBI.DAT'',/,/)')

 99   CONTINUE

      CLOSE(1,IOSTAT=IOS1)      

c      WRITE(*,FMT='(/,6X,''TRANSFERT MATRIX (Rij) IN THE COUPLED FRAME''
c     >,/)')
c      WRITE(*,100) ((R(I,J),J=1,4),I=1,4)
c 100     FORMAT(6X,4F13.8)


C     Computation of eigen values and eigen vectors
      
      REIG = R
      CALL RG(4,4,REIG,WR,WI,1,Q,GV1,FV1,IERR)         ! The output matrix Q is composed of as follows:
                                                       ! 1st column: real part of eigenvector 1
      EV1 = COMPLEX(WR(1),WI(1))                       ! 2nd column: imaginary part of eigenvector 1
      EV2 = COMPLEX(WR(2),WI(2))                       ! 3rd column: real part of eigenvector 2
      EV3 = COMPLEX(WR(3),WI(3))                       ! 4th column: imaginary part of eigenvector 2
      EV4 = COMPLEX(WR(4),WI(4))


C     Normalization of the eigen vectors


      NORM1  = SQRT(Q(1,1)**2+Q(1,2)**2+Q(2,1)**2+Q(2,2)**2+Q(3,1)**2+Q(
     >3,2)**2+Q(4,1)**2+Q(4,2)**2)
      NORM2  = SQRT(Q(1,3)**2+Q(1,4)**2+Q(2,3)**2+Q(2,4)**2+Q(3,3)**2+Q(
     >3,4)**2+Q(4,3)**2+Q(4,4)**2)

      V(1,1) = COMPLEX(Q(1,1),Q(1,2))/NORM1
      V(2,1) = COMPLEX(Q(2,1),Q(2,2))/NORM1
      V(3,1) = COMPLEX(Q(3,1),Q(3,2))/NORM1
      V(4,1) = COMPLEX(Q(4,1),Q(4,2))/NORM1
      V(1,2) = COMPLEX(Q(1,1),-Q(1,2))/NORM1
      V(2,2) = COMPLEX(Q(2,1),-Q(2,2))/NORM1
      V(3,2) = COMPLEX(Q(3,1),-Q(3,2))/NORM1
      V(4,2) = COMPLEX(Q(4,1),-Q(4,2))/NORM1
      V(1,3) = COMPLEX(Q(1,3),Q(1,4))/NORM2
      V(2,3) = COMPLEX(Q(2,3),Q(2,4))/NORM2
      V(3,3) = COMPLEX(Q(3,3),Q(3,4))/NORM2
      V(4,3) = COMPLEX(Q(4,3),Q(4,4))/NORM2
      V(1,4) = COMPLEX(Q(1,3),-Q(1,4))/NORM2
      V(2,4) = COMPLEX(Q(2,3),-Q(2,4))/NORM2
      V(3,4) = COMPLEX(Q(3,3),-Q(3,4))/NORM2
      V(4,4) = COMPLEX(Q(4,3),-Q(4,4))/NORM2


C     Phase uncertainity

      THETA1 = ATAN2(AIMAG(V(1,1)),REAL(V(1,1)))
      THETA2 = ATAN2(AIMAG(V(3,3)),REAL(V(3,3)))
      PHASE1 = COMPLEX(COS(-PI/2.D0-THETA1),SIN(-PI/2.D0-THETA1))
      PHASE2 = COMPLEX(COS(-PI/2.D0-THETA2),SIN(-PI/2.D0-THETA2))

      V(1,1) = COMPLEX(Q(1,1),Q(1,2))/NORM1*PHASE1
      V(2,1) = COMPLEX(Q(2,1),Q(2,2))/NORM1*PHASE1
      V(3,1) = COMPLEX(Q(3,1),Q(3,2))/NORM1*PHASE1
      V(4,1) = COMPLEX(Q(4,1),Q(4,2))/NORM1*PHASE1
      V(1,2) = COMPLEX(Q(1,1),-Q(1,2))/NORM1*CONJG(PHASE1)
      V(2,2) = COMPLEX(Q(2,1),-Q(2,2))/NORM1*CONJG(PHASE1)
      V(3,2) = COMPLEX(Q(3,1),-Q(3,2))/NORM1*CONJG(PHASE1)
      V(4,2) = COMPLEX(Q(4,1),-Q(4,2))/NORM1*CONJG(PHASE1)
      V(1,3) = COMPLEX(Q(1,3),Q(1,4))/NORM2*PHASE2
      V(2,3) = COMPLEX(Q(2,3),Q(2,4))/NORM2*PHASE2
      V(3,3) = COMPLEX(Q(3,3),Q(3,4))/NORM2*PHASE2
      V(4,3) = COMPLEX(Q(4,3),Q(4,4))/NORM2*PHASE2
      V(1,4) = COMPLEX(Q(1,3),-Q(1,4))/NORM2*CONJG(PHASE2)
      V(2,4) = COMPLEX(Q(2,3),-Q(2,4))/NORM2*CONJG(PHASE2)
      V(3,4) = COMPLEX(Q(3,3),-Q(3,4))/NORM2*CONJG(PHASE2)
      V(4,4) = COMPLEX(Q(4,3),-Q(4,4))/NORM2*CONJG(PHASE2)

 
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

      P(1,1) = DBLE(COMPLEX(0.D0,1.D0/SQRT(2.D0))*(V(1,1)-V(1,2)))
      P(2,1) = DBLE(COMPLEX(0.D0,1.D0/SQRT(2.D0))*(V(2,1)-V(2,2)))
      P(3,1) = DBLE(COMPLEX(0.D0,1.D0/SQRT(2.D0))*(V(3,1)-V(3,2)))
      P(4,1) = DBLE(COMPLEX(0.D0,1.D0/SQRT(2.D0))*(V(4,1)-V(4,2)))
      P(1,2) = DBLE(COMPLEX(1.D0/SQRT(2.D0),0.D0)*(V(1,1)+V(1,2)))
      P(2,2) = DBLE(COMPLEX(1.D0/SQRT(2.D0),0.D0)*(V(2,1)+V(2,2)))
      P(3,2) = DBLE(COMPLEX(1.D0/SQRT(2.D0),0.D0)*(V(3,1)+V(3,2)))
      P(4,2) = DBLE(COMPLEX(1.D0/SQRT(2.D0),0.D0)*(V(4,1)+V(4,2)))
      P(1,3) = DBLE(COMPLEX(0.D0,1.D0/SQRT(2.D0))*(V(1,3)-V(1,4)))
      P(2,3) = DBLE(COMPLEX(0.D0,1.D0/SQRT(2.D0))*(V(2,3)-V(2,4)))
      P(3,3) = DBLE(COMPLEX(0.D0,1.D0/SQRT(2.D0))*(V(3,3)-V(3,4)))
      P(4,3) = DBLE(COMPLEX(0.D0,1.D0/SQRT(2.D0))*(V(4,3)-V(4,4)))
      P(1,4) = DBLE(COMPLEX(1.D0/SQRT(2.D0),0.D0)*(V(1,3)+V(1,4)))
      P(2,4) = DBLE(COMPLEX(1.D0/SQRT(2.D0),0.D0)*(V(2,3)+V(2,4)))
      P(3,4) = DBLE(COMPLEX(1.D0/SQRT(2.D0),0.D0)*(V(3,3)+V(3,4)))
      P(4,4) = DBLE(COMPLEX(1.D0/SQRT(2.D0),0.D0)*(V(4,3)+V(4,4)))

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
         
      CALL TWSS(ALPHA1,BETA1,GAMMA1,ALPHA2,BETA2,GAMMA2,rPARAM,C,P)
      
      
C     Hamiltonian pertubation parameters and  Computing of the unperturbed tunes from the coupling parameters
      
      CALL HAMILT(NU1,NU2,rPARAM,CMOINS,DELTA,NUX0,NUY0,DELTA2,CPLUS)


C     Display of the results

      OPEN(2, FILE='ETparam.res',STATUS='UNKNOWN',IOSTAT=IOS2)

      nres = 7
c      CALL POSITI(NRES,
c     >                 ARCLEN,LABEL,KEYWOR)

      CALL EXTRAC(2,ARCLEN,R,rPARAM,C,NU1,NU2,ALPHA1,ALPHA2,BETA1,BETA2,
     >GAMMA1,GAMMA2,CMOINS,CPLUS,DELTA,DELTA2,NUX0,NUY0,P)
      
      CLOSE(2,IOSTAT=IOS2)


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
      
C        N_PROP = .TRUE.                                                    ! Flag turned on by default,
C        N_PROP = .FALSE.                                                    ! Flag turned on by default,
                                                                         ! it allows the computation of 
C        IF(N_PROP .EQV. .TRUE.) CALL PROPAG(3,4,7,8,P,C,NU1,NU2,CMOINS)    ! the generalized Twiss' parameters 
        CALL PROPAG(3,4,nres,8,P,C,NU1,NU2,CMOINS,NUX0,NUY0)    ! the generalized Twiss' parameters 
                                                                         ! and coupling strenght along 
        CLOSE(8,IOSTAT=IOS8)                                               ! the ring
        CLOSE(7,IOSTAT=IOS7)      
        CLOSE(4,IOSTAT=IOS4)
        CLOSE(3,IOSTAT=IOS3)
      ENDIF

      END
