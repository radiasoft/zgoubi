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
C  Brookhaven National Laboratory               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE MCOBJ3(KTIR,CENTRE,IMI,IMA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER  KTIR(*)*(*)
      DIMENSION CENTRE(*)
C     ----------------------------------------------------
C     Sorting at random inside three 2-D ellipses or other
C     ----------------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST2/ ZERO, UN
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/FITEX/ DN0,C0,C1,C2,C3,DL
      CHARACTER  KAR(41)
      COMMON/KAR/ KAR
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IMAXD,IMAXT
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/UNITS/ UNIT(MXJ)
 
      PARAMETER (MXJ1=MXJ-1)
      DIMENSION ALP(MXJ1),BET(MXJ1),EPS(MXJ1),RMA(MXJ1),RMB(MXJ1)

      DATA IKAR / 0 /

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
      RMA(2)=A(NOEL,53)
      RMA(4)=A(NOEL,63)
      RMA(6)=A(NOEL,73)
      RMB(2)=A(NOEL,54)
      RMB(4)=A(NOEL,64)
      RMB(6)=A(NOEL,74)
 
C----- LECTURE GENERATEURS ( SEULEMENT AU 1-ER PASSAGE SI REBELOTE)
C      TIRAGE AVEC MELANGE ALEATOIRE
      IF(IPASS .EQ. 1) THEN
        IR1 = A(NOEL,80)
        IR2 = A(NOEL,81)
        IR3 = A(NOEL,82)
        IR1=(IR1/2)*2+1
        IR2=(IR2/2)*2+1
        IR3 =(IR3 /2)*2+1
      ENDIF
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,100)
100     FORMAT(/,15X,' Distribution  with  ellipcal  frontiers')
 
        WRITE(NRES,123) (CENTRE(J),J=2,6),CENTRE(1)
 123    FORMAT(/,15X,' Ellipse  centers  (MKSA units) : '
     >  ,/,11X,' horizontal    ( Yo, To ):',T50,1P,2G12.4
     >  ,/,11X,' vertical      ( Zo, Po ):',T50,   2G12.4
     >  ,/,11X,' longitudinal  ( Xo, Do ):',T50,   2G12.4,/)
 
        WRITE(NRES,109) (ALP(J),BET(J),EPS(J),RMA(J),RMB(J),J=2,6,2)
 109    FORMAT(/,15X,
     >  ' Alpha (rad), Beta(m/rad), E/pi(m.rad), Cut-off(*BetEps/pi) :'
     >  ,/,11X,' horizontal    :',T35,1P,3G12.4,G10.2,' (',G10.2,')'
     >  ,/,11X,' vertical      :',T35,   3G12.4,G10.2,' (',G10.2,')'
     >  ,/,11X,' longitudinal  :',T35,   3G12.4,G10.2,' (',G10.2,')')

        WRITE(NRES,FMT='(/,15X,'' Sorting  types  (Y/Z/L) : '',
     >      3(A9,''/''))') (KTIR(J),J=2,MXJ1,2)
      ENDIF
   
C----- Constitution of the beam
C MXJ1=6 ; J : 2,4,6 -> Y,Z,s 
      DO 1 J=2,MXJ1,2
C------- J1 : 3,5,7(1) -> T,P,D
        J1=J+1
        IF(J1.EQ.MXJ) J1=1 
        IF(EPS(J).EQ.0.D0) THEN
          DO I=IMI,IMA
            FO(J ,I)=0.D0
            FO(J1,I)=0.D0
          ENDDO
        ELSE
          REB=SQRT(EPS(J)*BET(J))
          REBM=SQRT(RMA(J)*EPS(J)*BET(J))
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
     >            SQRT(RMA(J)*EPS(J)*(1.D0+ALP(J))**2/BET(J))
                SIGN =  1.D0 
                IF(2.D0*(RNDM()-.5D0).LE.0.D0) SIGN=-SIGN
                X = R*SIGN
                FO(J1,I) = X/UNIT(J)
              ENDIF

            ELSEIF(KTIR(J) .EQ. 'Gaussian') THEN
C-------------  Tirage uniforme en exp(-r2) = gaussien en y et y'

              SM=EXP(-RMA(J)*RMA(J)/2.D0)
C              SM=EXP(-RMA(J)/2.D0)
              IF(RMB(J) .EQ. 0.D0) THEN
C---------------- Sorting in [0,SMA]
C                R=RNDM()*(1.D0-SM)+SM
                R=1.D0 + RNDM()*(SM-1.D0)
              ELSE
C---------------    Sorting in [SMA,SMB]
CCCCC                SMB = EXP(-RMB(J)/2.D0)
                SMB = EXP(-RMB(J)*RMB(J)/2.D0)
                SM  = EXP(-RMA(J)*RMA(J)/2.D0)
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
              DXM = AXB*AXB - (GAM*X*X - EPS(J))/BET(J)
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
                  
      DO 4 I=IMI,IMA
        FO(1,I)=FO(1,I) + CENTRE(1)/UNIT(6)
        DO 4 J=2,MXJ1
          FO(J,I)=FO(J,I) + CENTRE(J)/UNIT(J-1)
 4      CONTINUE
 
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
      END
