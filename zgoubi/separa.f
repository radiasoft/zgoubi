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
C  Upton, NY, 11973
C  -------
      SUBROUTINE SEPARA(*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH

      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      COMMON/DON/ A(09876,99),IQ(09876),IP(09876),NB,NOEL
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.MARK.H"     ! COMMON/MARK/ KART,KALC,KERK,KUASEX
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AMASS,Q,G,TOO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE 'MXFS.H'
      INCLUDE 'MXSCL.H'
      INCLUDE "C.SCAL.H"     ! COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
C      COMMON/SCAL/SCL(MXF,MXS),TIM(MXF,MXS),NTIM(MXF),JPA(MXF,MXP),KSCL
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)

      CHARACTER(10) TYP(2)
      DATA TYP / 'HORIZONTAL',' VERTICAL ' /

C     ... FACTEUR D'ECHELLE DES ChampS. UTILISE PAR 'SCALING'
      SCAL = SCAL0()
      IF(KSCL .EQ. 1) SCAL = SCAL0()*SCALER(IPASS,NOEL,
     >                                                 ZERO)
C     ...........................................................
C        IJ = 0 ( DRIFT SPACE) , 1 (H SEPARATION) OR 2 ( V SEP.)
C        DL = LONGUEUR ( M )
C        E  ( V/M ) = E-FIELD  (USUALLY > 0)
C        B  ( T )   = B-FIELD  (USUALLY < 0)
C     ...........................................................

      IJ = NINT( A(NOEL,1) )
      DL = A(NOEL,2)
      E =  A(NOEL,3)
      B =  A(NOEL,4)

      XL = DL*1.D2

      IF(IJ .NE. 0) THEN
C        P0 = BORO*CL*1.D-9 *Q/QE
C        P0 = BORO*CL*1.D-9 *Q
        P0 = BORO*CL9 *Q
        IF(NRES.GT.0) THEN
          WRITE(NRES,110)
     >      TYP(IJ),DL,E,B,AMASS,Q,P0/SQRT(P0*P0+AMASS*AMASS)
 110      FORMAT(/,25X,' ******  ELECTROSTATIC  MASS  SEPARATOR  ******'
     >          ,/,40X,A
     >          ,/,30X,' Length                    = ',1P,G12.5,'  m'
     >          ,/,30X,' E                         = ',G12.5,' V/m '
     >          ,/,30X,' B                         = ',G12.5,'  T  '
     >        ,/,/,30X,' Mass  of  particles       = ',G12.5,' MeV/c2'
     >        ,/,/,30X,' Charge  of  particles     = ',G12.5,' C'
     >          ,/,30X,' Reference  beta           = ',G12.5)
          IF(B .NE. 0.D0) WRITE(NRES,111) E/(-B*CL)
 111      FORMAT(30X,' Beta  wanted = -E/v.B =',1P,G12.5,//)
        ENDIF

C        IF(AMASS*Q .EQ. 0.D0) THEN
        IF(Q*AMASS.EQ.0D0) THEN
          WRITE(NRES,106)
 106      FORMAT(//,15X,' Please  give   M  &  Q  of  projectiles !'
     >           ,/,15X,' - use  keyword  ''PARTICUL''',/)
          RETURN 1
        ENDIF

C        IF(ABS(Q/QE) .GE. 2.D0) THEN
        IF(ABS(Q) .GE. 2.D0) THEN
          WRITE(NRES,FMT='(15X,'' SORRY, SUBROUTINE SEPARA '',
     >    ''WORKS ONLY FOR SINGLE-CHARGE PARTICLES'',/)')
          RETURN 1
        ENDIF
      ELSEIF(IJ .EQ. 0.D0) THEN

        IF(NRES.GT.0) WRITE(NRES,109) DL
 109    FORMAT(//,25X,' +++++++++++ SEPARATOR OFF ++++++++++++'
     >  ,//,22X,' Equivalent  to  drift  with  length ',F12.5,' m',//)
C     >  ,//,22X,' EQUIVALENT A UN ESPACE  LIBRE  DE',F12.5,' M',//)

      ENDIF

C----- OPTION COMPTAGE DES TRAJECTOIRES
C      HORS DES LIMITES   Y - Z  DU  SEPARATEUR
      IF(LIMIT .EQ. 1) THEN
        IF(XL.GE.0.D0) CALL CHMBR(1,IMAX)
      ENDIF

      IF(IJ .NE. 0) THEN

        DO 10 I=1,IMAX

          IF    (IJ .EQ. 1) THEN
C           ** HORIZONTAL  SEPARATION
            XO = F(2,I)*.01D0
            TO = F(3,I)*.001D0
            ZO = F(4,I)*.01D0
            PO = F(5,I)*.001D0
          ELSEIF(IJ .EQ. 2) THEN
C           ** VERTICAL  SEPARATION
            ZO = F(2,I)*.01D0
            PO = F(3,I)*.001D0
            XO = F(4,I)*.01D0
            TO = F(5,I)*.001D0
          ENDIF

          P = P0*F(1,I)
          BTA = P/SQRT(P*P+AMASS*AMASS)
          G1 = SQRT(1-BTA*BTA)
          V = BTA*CL


          C1 = V*SIN(TO)*COS(PO)
          C2 = V*COS(TO)*COS(PO)
          ZP= V*SIN(PO)
          IF( E*B .NE. 0.D0) THEN
            ALP = -E*CL*CL/(AMASS*1.D6)*G1
            OME = -B*CL*CL/(AMASS*1.D6)*G1
            R = SQRT( C1*C1 + (C2+ALP/OME)**2 ) / OME
            IF(C1 .EQ. 0.D0) C1=1.D-10
            EPS = ATAN2( C2+ALP/OME , C1 )
            CALL ZEROF(I,DL,BTA,ALP,OME,C1,R,EPS,T)
            ARG = OME*T+EPS
            RCA = R*COS(ARG)
            RSA = R*SIN(ARG)
            X =  RSA - (ALP/OME+C2)/OME + XO
            Y = -RCA - (ALP*T  -C1)/OME
            Z = ZP*T + ZO
            XP = OME*RCA
            YP = OME*RSA - ALP/OME
            THET = ATAN2(XP,YP)
            PHI  = ATAN2(ZP,SQRT(XP*XP+YP*YP))
          ELSEIF( B .EQ. 0.D0) THEN
            ALP = -E*CL*CL/(AMASS*1.D6)*G1
            T = DL/V
            X =  (-ALP*T*.5D0 + C1)*T + XO
            XP=  -ALP*T + C1
            YP = C2
            Z = ZP*T + ZO
            THET = ATAN2(XP,YP)
            PHI  = ATAN2(ZP,SQRT(XP*XP+YP*YP))
          ELSEIF( E .EQ. 0.D0) THEN
            OME = -B*CL*CL/(AMASS*1.D6)*G1
            R = SQRT( C1*C1 + C2*C2 ) / OME
            IF(C1 .EQ. 0.D0) C1=1.D-10
            EPS = ATAN2( C2 , C1 )
            T = ( ACOS( (C1/OME - DL)/R ) - EPS )/OME
            ARG = OME*T+EPS
            RCA = R*COS(ARG)
            RSA = R*SIN(ARG)
            X =  RSA - C2/OME + XO
            Y = -RCA + C1/OME
            Z = ZP*T + ZO
            XP = OME*RCA
            YP = OME*RSA
            THET = ATAN2(XP,YP)
            PHI  = ATAN2(ZP,SQRT(XP*XP+YP*YP))
          ENDIF

          IF    (IJ .EQ. 1) THEN
C           ** HORIZONTAL  SEPARATION
            F(2,I) = X*1.D2
            F(3,I) = THET*1.D3
            F(4,I) = Z*1.D2
            F(5,I) = PHI*1.D3
          ELSEIF(IJ .EQ. 2) THEN
C           ** VERTICAL  SEPARATION
            F(2,I) = Z*1.D2
            F(3,I) = PHI*1.D3
            F(4,I) = X*1.D2
            F(5,I) = THET*1.D3
          ENDIF
          F(6,I) = F(6,I) + V*T*1.D2
          F(7,I) = F(7,I) + T   *1.D+6

 10     CONTINUE

C----- CASSURE DE Z-SYMETRIE EVENTUELLE
        ZSYM=.FALSE.

        IF( KSPN .EQ. 1 ) THEN
          IF    (IJ .EQ. 1) THEN
            DO I=1,IMAX
              SX = SF(1,I)
              SY = SF(2,I)
              CALL ROTZ(-G/G1*OME*T,SX,SY)
              SF(1,I) = SX
              SF(2,I) = SY
            ENDDO
          ENDIF
        ENDIF

      ELSEIF(IJ .EQ. 0) THEN
C------- SEPARATEUR OFF
        CALL ESL(0,XL,1,IMAX)

      ENDIF

C----- OPTION COMPTAGE DES TRAJECTOIRES
C      HORS DES LIMITES   Y - Z  DU  SEPARATEUR
      IF(LIMIT .EQ. 1) THEN
        IF(XL.GE.0.D0) CALL CHMBR(1,IMAX)
      ENDIF

      RETURN
      END
