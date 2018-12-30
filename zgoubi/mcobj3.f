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
C  François Méot <fmeot@bnl.gov>
C  Brookhaven National Laboratory
C  C-AD, Bldg 911
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE MCOBJ3(KTIR,CENTRE,KNRM,IMI,IMA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*)  KTIR(*)
      DIMENSION CENTRE(*)
C     ----------------------------------------------------
C     Sorting at random inside three 2-D ellipses or other
C     ----------------------------------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
      CHARACTER(1) LET
      INCLUDE "C.FAISCT.H"     ! COMMON/FAISCT/ LET(MXT)
      INCLUDE "C.FITEX.H"     ! COMMON/FITEX/ DN0,C0,C1,C2,C3,DL
      CHARACTER(LEN=1) KAR(41)
      INCLUDE "C.KAR.H"     ! COMMON/KAR/ KAR
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)

      PARAMETER (MXJ1=MXJ-1)
      DIMENSION ALP(MXJ1),BET(MXJ1),EPS(MXJ1),RMA(MXJ1),RMB(MXJ1)
      DIMENSION EPST(MXJ1)

      DIMENSION DIS(MXJ)
      SAVE IRANK

      INTEGER DATE_TIME (8)
      CHARACTER (LEN = 12) REAL_CLOCK (3)

      DATA IKAR / 0 /
      DATA IRANK / -1 /

C     .. PARAMETRE ELLIPSES
      ALP(2)=A(NOEL,50)
      ALP(4)=A(NOEL,60)
      ALP(6)=A(NOEL,70)
      BET(2)=A(NOEL,51)
      BET(4)=A(NOEL,61)
      BET(6)=A(NOEL,71)
      EPS(2)=A(NOEL,52)
      EPS(4)=A(NOEL,62)
      EPS(6)=A(NOEL,72)
C      RMA(2)=A(NOEL,53)
C      RMA(4)=A(NOEL,63)
C      RMA(6)=A(NOEL,73)
      J = 1
      DOWHILE (J .LE. 3)
        RMA(2*J)=A(NOEL,53 + 10*(J-1))
        IF(RMA(2*J) .LT. 0.D0) THEN
          RMA(2*J) = -RMA(2*J)
          RMB(2*J) = A(NOEL,54 + 10*(J-1))
          DIS(2*J)  = A(NOEL,55 + 10*(J-1))
          DIS(2*J+1) = A(NOEL,56 + 10*(J-1))
        ELSE
          DIS(2*J)  = A(NOEL,54 + 10*(J-1))
          DIS(2*J+1) = A(NOEL,55 + 10*(J-1))
        ENDIF
C      RMB(4)=A(NOEL,64)
C      RMB(6)=A(NOEL,74)
        J = J + 1
      ENDDO

      IF    (KNRM .EQ. 1) THEN
        EPST(2)=EPS(2) / CENTRE(1)
        EPST(4)=EPS(4) / CENTRE(1)
        EPST(6)=EPS(6) / CENTRE(1)
      ELSE
        EPST(2)=EPS(2)
        EPST(4)=EPS(4)
        EPST(6)=EPS(6)
      ENDIF

C----- LECTURE GENERATEURS ( SEULEMENT AU 1-ER PASSAGE SI REBELOTE)
C      TIRAGE AVEC MELANGE ALEATOIRE
      IF(IPASS .EQ. 1) THEN
C        IR1 = NINT(A(NOEL,80))
C        IR2 = NINT(A(NOEL,81))
C        IR3 = NINT(A(NOEL,82))
C FM, VR. 2017
C ModifIED FOR MULTIPLE-RUN ON NERSC
        IF(IRANK .EQ. -1) THEN
          IR1 = 2*(NINT(A(NOEL,80))/2)+1
        ELSE
C GENERATE A RANDOM SEED TO INITIAL THE SERIES
C FIRST, INITIALIZE A SERIES
          CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2),
     >                 REAL_CLOCK (3), DATE_TIME)
Consider the following example executed on 2000 March 28 at 11:04:14.5:
c     This assigns the value "20000328" to REAL_CLOCK (1), the value "110414.500" to REAL_CLOCK (2),
C     and the value "-0500" to REAL_CLOCK (3). The following values are assigned to
C     DATE_TIME: 2000, 3, 28, -300, 11, 4, 14, and 500.
          READ(REAL_CLOCK (2)(5:12),*) SSI
          IR1 = NINT(SSI*1000.D0)
        ENDIF
        TEMP = RNDM2(IR1)
      ENDIF

      IF(NRES.GT.0) THEN
        WRITE(NRES,100)
100     FORMAT(/,15X,' Distribution  with  ellipcal  frontiers')

        WRITE(NRES,123) (CENTRE(J),J=2,6),CENTRE(1)
 123    FORMAT(/,15X,' Ellipse  centers  (MKSA units) : '
     >  ,/,11X,' horizontal    ( Yo, To ):',T50,1P,2G12.4
     >  ,/,11X,' vertical      ( Zo, Po ):',T50,   2G12.4
     >  ,/,11X,' longitudinal  ( Xo, Do ):',T50,   2G12.4,/)

        IF(KNRM .EQ. 1)
     >  WRITE(NRES,FMT='(A,/)') ' Transverse emittances as declared are'
     >  //' assumed normalized to BORO. Actual geometrical emittances '
     >  //' below drawn by changing  epsilon -> epsilon/(central p/p0).'

        WRITE(NRES,109) (ALP(J),BET(J),EPST(J),RMA(J),RMB(J),J=2,6,2)
 109    FORMAT(/,15X,
     >  ' Alpha (rad), Beta(m/rad), E/pi(m.rad), Cut-off(*BetEps/pi) :'
     >  ,/,11X,' horizontal    :',T35,1P,3G12.4,G10.2,' (',G10.2,')'
     >  ,/,11X,' vertical      :',T35,   3G12.4,G10.2,' (',G10.2,')'
     >  ,/,11X,' longitudinal  :',T35,   3G12.4,G10.2,' (',G10.2,')')

        WRITE(NRES,108) (DIS(J),J=2,5)
 108    FORMAT(/,15X
     >  ,' Coordinates correlated :  '
     >  ,' Y ->  Y +D_Y * dp/p,  T ->  T +D_T * dp/p,  '
     >  ,' Z -> ZY +D_Z * dp/p,  P ->  P +D_P * dp/p.  '
     >  ,/,11X,' Dispersion values as follows :'
     >  ,/,11X,' horizontal, D_Y, D_T    :',T35,1P,2(G12.4,2X)
     >  ,/,11X,' vertical, D_Z and D_P   :',T35,1P,2(G12.4,2X))

        WRITE(NRES,FMT='(/,15X,'' Sorting  types  (Y/Z/L) : '',
     >      3(A9,''/''))') (KTIR(J),J=2,MXJ1,2)

        IF(IRANK.NE.-1) WRITE(NRES,FMT='(/,A,I0)')
     >  'IR!-2 updated using IRANK =  "',IRANK,'"'

      ENDIF

C----- Constitution of the beam
C MXJ1=6 ; J : 2,4,6 -> Y,Z,s
      DO 1 J=2,MXJ1,2
C------- J1 : 3,5,7(1) -> T,P,D
        J1=J+1
        IF(J1.EQ.MXJ) J1=1
        IF(EPST(J).EQ.0.D0) THEN
          DO I=IMI,IMA
            FO(J ,I)=0.D0
            FO(J1,I)=0.D0
          ENDDO
        ELSE
          REB=SQRT( EPST(J)*BET(J) )
          REBM=SQRT(RMA(J)*EPST(J)*BET(J))
C          RM = RMA(J)*REB
          DO I=IMI,IMA

            IF    (KTIR(J) .EQ. 'Uniform') THEN
C       Tirage uniforme en r2 = uniforme en surface, elliptique en y et y'

              IF(J.LE.4) THEN
C                TRANSVERSE COORDINATES
                R=SQRT(RNDM())*REBM
                ANG = 2.D0*(RNDM()-.5D0)*PI
                X = R*COS(ANG)
                FO(J ,I) = X/UNIT(J-1)
                FO(J1,I) = (R*SIN(ANG)-ALP(J)*X)/BET(J)/UNIT(J)
              ELSE
C                DP/P, X
                R=RNDM()*REBM
                SIGN =  1.D0
                IF(2.D0*(RNDM()-.5D0).LE.0.D0) SIGN=-SIGN
                X = R*SIGN
                FO(J ,I) = X/UNIT(J-1)
                R=RNDM()*
     >            SQRT(RMA(J)*EPST(J)*(1.D0+ALP(J))**2/BET(J))
                SIGN =  1.D0
                IF(2.D0*(RNDM()-.5D0).LE.0.D0) SIGN=-SIGN
                X = R*SIGN
                FO(J1,I) = X/UNIT(J)
              ENDIF

            ELSEIF(KTIR(J) .EQ. 'Gaussian') THEN
C-------------  Tirage uniforme en exp(-r2) = gaussien en y et y'

              SM = EXP(-RMA(J)*RMA(J)/2.D0)
              IF(RMB(J) .EQ. 0.D0) THEN
C---------------- Sorting in [0,SMA]
                R=1.D0 + RNDM()*(SM-1.D0)
              ELSE
C---------------    Sorting in [SMA,SMB]
                SMB = EXP(-RMB(J)*RMB(J)/2.D0)
                R=SM + RNDM()*(SMB-SM)
              ENDIF

CCCCC------------- /2.D0 leads to same emittance as calculated from
CCCCC                            concentration ellipse in zgplot :
CCCCC             R=SQRT(-2.D0*LOG(R))*REB /2.D0

C-------------- Leads to sigma_x=sqrt( beta_x epsilon_x/pi) :
              R=SQRT(-2.D0*LOG(R))*REB
              ANG = 2.D0*(RNDM()-.5D0)*PI
              X = R*COS(ANG)
              FO(J ,I) = X/UNIT(J-1)
              FO(J1,I) = (R*SIN(ANG)-ALP(J)*X)/BET(J)/UNIT(J)

            ELSEIF(KTIR(J) .EQ. 'Parabolic') THEN
C       Tirage parabolic en y : p(y) = (1-y2/y02)*3/4/y0

 7            R = RNDM()
              A3 = ACOS( 1.D0 - 2.D0 * R ) / 3.D0
              X = - REBM * COS ( A3 + PI/3.D0 )
              FO(J ,I) = X/UNIT(J-1)

              AXB = ALP(J)*X/BET(J)
              GAM = (1.D0 + ALP(J)*ALP(J))/ BET(J)
              DXM = AXB*AXB - (GAM*X*X - EPST(J))/BET(J)
              IF(DXM .LT. 0.D0) GOTO 7
              DXM = SQRT( DXM )
              XM = -AXB
              R = RNDM()
              A3 = ACOS( 1.D0 - 2.D0 * R ) / 3.D0
              X = XM - DXM * COS ( A3 + PI/3.D0 )
              FO(J1,I) = X/UNIT(J)
            ENDIF
          ENDDO
        ENDIF
 1    CONTINUE

      DO I=IMI,IMA
        DO J=2,MXJ1-1
c       write(88,fmt='(a,1p,2(i4,2x),6(e12.4,2x))') ' mcobj3 ',
c     > i,j,FO(J,I),DIS(J) ,FO(1,I),DIS(J)/UNIT(J-1) * FO(1,I),
c     >  FO(J,I) + CENTRE(J)/UNIT(J-1),
c     >   FO(J,I) + CENTRE(J)/UNIT(J-1) + DIS(J)/UNIT(J-1) * FO(1,I)
c              read(*,*)
          FO(J,I)=FO(J,I) + DIS(J)/UNIT(J-1) * FO(1,I)
          FO(J,I)=FO(J,I) + CENTRE(J)/UNIT(J-1)
        ENDDO
        J = 1
        FO(J,I)=FO(J,I) + CENTRE(J)/UNIT(6)
C        FO(1,I)=FO(1,I) + CENTRE(1)/UNIT(6)
        J = MXJ1
        FO(J,I)=FO(J,I) + CENTRE(J)/UNIT(J-1)
      ENDDO

      ENTRY MCOBJA

      DO I=1,IMAX
        F(1,I)=FO(1,I)
        DO J=2,MXJ
          F(J,I)=FO(J,I)
        ENDDO
        IKAR=IKAR+1
        IF(IKAR.GT.41) IKAR=1
        LET(I)= KAR(IKAR)
        IREP(I)=I
        IEX(I)=1
      ENDDO

      RETURN

      ENTRY MCOBJB(IRKI)
      IRANK = IRKI
      RETURN

      END
