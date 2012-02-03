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
      SUBROUTINE AIMANT(ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------
C     Optical elements defined in polar frame
C-----------------------------------------------------
      COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
C----- Conversion  coord. (cm,mrd) -> (m,rd)
      INCLUDE 'MAXCOO.H'
      COMMON/UNITS/ UNIT(MXJ)
 
      LOGICAL BONUL

      PARAMETER (I0=0, ZERO=0.D0)

      PARAMETER (MSR=10)
      CHARACTER QSHROE(MSR)*(2),QSHROS(MSR)*(2)
      DIMENSION VSHROE(MSR),    VSHROS(MSR)

C------ Fields defined in polar coordinates - Champs definis en coordonnees polaires
      KART = 2
 
      CALL CHXP(ND,KALC,KUASEX,
     >                         XL,DSREF,NDD)
      IF(NRES .GT. 0) CALL FLUSH2(NRES,.FALSE.)

      IF ( BONUL(XL,PAS) ) GOTO 99

      CALL SCUMW(DSREF)
      CALL SCUMR(
     >            DUM,SCUM,TCUM) 

      CALL RAZDRV(3)

      IF(PAS .LT. 0.D0)  THEN
        TEMP=AE
        AE=-AS
        AS=-TEMP
        TEMP=XI
        XI=XF
        XF=TEMP
      ENDIF

C----- DXI = step angle in polar frame
      DXI=PAS/RM

C      IF(KP .NE. 2)  THEN
      IF(KP .EQ. 1)  THEN
        DP = A(NOEL,ND+20)
        IF(DP.LE.0.D0)  CALL INITRA(1)
        IF(NRES.GT.0) WRITE(NRES,100) DP
        X=AE+AT/2.D0
        Y=RM
        Z=ZERO
        T=ZERO
        P=ZERO
        XLIM=XI
        DXI=-DXI
        PAS=-PAS
        CALL INTEGR(.FALSE.,.TRUE.,.FALSE.)
        AA=-AE
        CALL CHAREF(.FALSE.,ZERO,ZERO,AA)
        RE=Y
        TE=T
        X=AE+AT/2.D0
        Y=RM
        T=ZERO
        XLIM=XF
        DXI=-DXI
        PAS=-PAS
        CALL INTEGR(.FALSE.,.TRUE.,.FALSE.)
        CALL CHAREF(.FALSE.,ZERO,ZERO,AS)
        RS=Y
        TS=T
C        GOTO 3
C      ENDIF
      ELSEIF(KP .EQ. 2) THEN
        RE = A(NOEL,NDD)
        TE = A(NOEL,NDD+1)
        RS = A(NOEL,NDD+2)
        TS = A(NOEL,NDD+3)

        IF(PAS .LT. 0.D0) THEN
          TEMP=TE
          TE=TS
          TS=TEMP
          TEMP=RE
          RE=RS
          RS=TEMP
        ENDIF
      ENDIF

C3     CONTINUE
      IF(NRES.GT.0) WRITE(NRES,101)RE,TE,RS,TS

      ALE=AE-TE
      XCE=-RE*SIN(TE)
      YCE=-RE*COS(TE)
      ALS=AS+TS
      XCS=-RS*SIN(AS)
      YCS=RS*COS(AS)
      X=XI
      XLIM=XF
      KP=2

      qshroE(1) = 'XS'
      VSHROE(1) = XCE
      qshroE(2) = 'YS'
      VSHROE(2) = YCE
      qshroE(3) = 'ZR'
      VSHROE(3) = ALE
      VSHROE(MSR) = 3
      qshroS(1) = 'XS'
      VSHROS(1) = XCS
      qshroS(2) = 'YS'
      VSHROS(2) = YCS
      qshroS(3) = 'ZR'
      VSHROS(3) = ALS
      VSHROS(MSR) = 3

      CALL TRANSF(QSHROE,VSHROE,QSHROS,VSHROS)

 99   CONTINUE
C----- Unset wedge correction, in case it has been set by MULTIPOL, BEND, etc.
      CALL INTEG3
C----- Unset coded step
      CALL CHXC1W(I0,I0)
      CALL DEPLAW(.FALSE.,I0)

      IF(NRES .GT. 0)
     >WRITE(NRES,FMT='(/,'' Cumulative length of optical axis = '',
     >1P,G17.9,
     >'' m ;  Time  (for ref. rigidity & particle) = '', 
     >1P,G14.6,'' s '')')  SCUM*UNIT(5), TCUM

      RETURN
 
100   FORMAT(/,5X,'  Magnet positionning computed from'
     > ,1X,'backward tracking of particle with impulse D =',F7.4,/)
101   FORMAT(/,30X,' Position of reference orbit on mechanical ',
     > ' faces',/,40X,' at entrance    RE =',1P,G14.6,' cm  TE =',
     > G14.6,' rad',/,40X,
     >     ' at exit        RS =',G14.6,' cm  TS =',G14.6,' rad',/)
C100   FORMAT(/,5X,'  Positionnement de l''aimant calculee par retour'
C     > ,1X,'arriere de la trajectoire d''impulsion D =',F7.4,/)
C101   FORMAT(/,30X,' Position de l''orbite moyenne sur les faces',
C     > ' mecaniques',/,40X,' d''ENTREE  RE =',F8.3,' cm  TE =',
C     > F9.6,' rad',/,40X,
C     > 'de SORTIE  RS =',F8.3,' cm  TS =',F9.6,' rad',/)
      END
