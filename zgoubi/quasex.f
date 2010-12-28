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
C  USA
C  -------
      SUBROUTINE QUASEX(ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------
C     Optical elements defined in cartesian frame
C-----------------------------------------------------
      COMMON/AIM/ BO, RO,FG,  GF,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
C      COMMON/CHAFUI/ XE,XS,CE(6),CS(6),QCE(6),QCS(6)
C FM, Oct 2001      
Ccccc COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/CONST/ CL9,CEL,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      COMMON/TRNSF/ XFE,XFS
C----- Conversion  coord. (cm,mrd) -> (m,rd)
      INCLUDE 'MAXCOO.H'
      COMMON/UNITS/ UNIT(MXJ)
 
      LOGICAL BONUL, LNUL
      LOGICAL MIRROR, LDUM

      PARAMETER(I0=0, ZERO=0.D0, IZERO=0)

C FIELDS ARE DEFINED IN CARTESIAN COORDINATES
      KART = 1
      CALL CHXC(ND,KALC,KUASEX,BORO,
     >                              XL,DSREF)
      IF(NRES .GT. 0) CALL FLUSH2(NRES,.FALSE.)
 
      IF ( LNUL(XL) ) RETURN
      IF ( BONUL(XL,PAS) ) RETURN

      CALL SCUMW(DSREF)
      CALL SCUMR(
     >           DUM,SCUM,TCUM) 

      CALL RAZDRV(3)

      DXI=PAS
      XFE=-XE
      XFS=XS-XLIM
      IF(KP .EQ. 1)  THEN
C------- Optical elmnt undergoes new positionning. Its frame becomes the new 
C        reference frame. No exit CHANGREF needed
        XCS = 0.D0
        YCS = 0.D0
        ALS = 0.D0
      ELSEIF(KP .EQ. 2)  THEN
C------- KP=2, element is misaligned. Hence compute XCS,YCS,ALS for automatic 
C        re-positionning of the exit frame to old position:
C        XLIM = XF
C        CL=COS(ALE)
C        SL=SIN(ALE)
C        XTEMP=XCE-XLIM*(1.D0-CL)
C        YTEMP=YCE+XLIM*SL
C Modified Jan 06, FM
C         WRITE(*,*) ' QUASEX ',XE,XS,XLIM
        XLM = XS-XE
        CL=COS(ALE)
        SL=SIN(ALE)
        XTEMP=XCE-XLM*(1.D0-CL)
        YTEMP=YCE+XLM*SL

        IF(PAS.GT.0.D0)  GOTO 1
        TEMP=XCE
        XCE=XTEMP
        XTEMP=TEMP
        TEMP=YCE
        YCE=YTEMP
        YTEMP=TEMP
 1      CONTINUE
        XCS=-XTEMP*CL-YTEMP*SL
        YCS=XTEMP*SL-YTEMP*CL
        ALS=-ALE
        IF(NRES.GT.0) WRITE(NRES,100) XCE,YCE,ALE
      ELSEIF(KP .EQ. 3)  THEN
C------- Optical elmnt is Z-rotated. Entrance and exit frames are 
C        tilted by (normally) half the deviation. 
C Modified, FM, Dec 05
C        YCS = -YCE
C        YCE = 0.D0
        XCS = 0.D0
        YCS = -YCE/COS(ALE)
        ALS=ALE
        CALL TRANSR(MIRROR,LDUM)
        IF(MIRROR) THEN 
          ALS = -(PI - ALE)
        ENDIF
      ENDIF

      IF(PAS.GT.0.D0)  GOTO 3
      TEMP=XI
      XI=XF
      XF=TEMP
      XFE=XLIM-XS
      XFS=XE
3     CONTINUE
      X=XI
      XLIM=XF
 
      CALL TRANSF

C----- Unset wedge correction, in case it has been set by MULTIPOL, BEND, etc.
      CALL INTEG3

      IF(NRES .GT. 0)
     >WRITE(NRES,FMT='(/,'' Cumulative length of optical axis = '',
     >1P,G17.9,
     >'' m ;   corresponding Time  (for ref. rigidity & particle) = '', 
     >1P,G14.6,'' s '')')  SCUM*UNIT(5), TCUM

      RETURN
C100   FORMAT(/,5X,'ELEMENT  DECENTRE  PAR  RAPPORT  A  L''AXE  OPTIQUE'
C     1,/,10X,'CENTRE  DE  LA  FACE  D''ENTREE  EN  X =',1P,G10.2
C     2,' CM   Y =',G10.2,' CM   INCLINAISON =',G12.4,' RAD',/)
100   FORMAT(/,5X,'Element  is  mis-aligned  wrt.  the  optical  axis'
     1,/,10X,'Center  of  entrance  EFB  is  at    X =',1P,G12.4
     2,' CM   Y =',G12.4,' cm,  tilt  angle =',G14.6,' RAD',/)
      END
