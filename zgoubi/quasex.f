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
      SUBROUTINE QUASEX(ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------
C     Optical elements defined in cartesian frame
C-----------------------------------------------------
      INCLUDE "C.AIM.H"     ! COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
                            ! COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      INCLUDE "C.CHAFUI.H"     ! COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
      INCLUDE "C.CONST_3.H"      ! COMMON/CONST/ CL9,CEL,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      INCLUDE "C.MARK.H"     ! COMMON/MARK/ KART,KALC,KERK,KUASEX
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE "C.SYNRA.H"     ! COMMON/SYNRA/ KSYN
      INCLUDE "C.TRAJ.H"     ! COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      INCLUDE "C.TRNSF.H"     ! COMMON/TRNSF/ XFE,XFS
C----- Conversion  coord. (cm,mrd) -> (m,rd)
      INCLUDE 'MAXCOO.H'
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)
 
      LOGICAL BONUL, LNUL
      LOGICAL MIRROR, LDUM

      PARAMETER(I0=0, ZERO=0.D0)
      SAVE DTTA, ZCE, PHI

      PARAMETER (MSR=8)
      CHARACTER(2) QSHROE(MSR),QSHROS(MSR)
      DIMENSION VSHROE(MSR), VSHROS(MSR)
      LOGICAL OKPRSR

      INTEGER DEBSTR, FINSTR
      CHARACTER(51) TXTF(4)
      DATA TXTF / 
     >'Change  of  frame  at  exit  of  element.',
     >'Element  is  mis-aligned  wrt.  the  optical  axis.',
     >'Automatic  positionning  of  element.',
     >'Special  alignement.'
     >/

C FIELDS ARE DEFINED IN CARTESIAN COORDINATES
      KART = 1
      CALL CHXC(ND,KALC,KUASEX,BORO,DPREF,
     >                              XL,DSREF,QSHROE,VSHROE)
      IF(NRES .GT. 0) CALL FLUSH2(NRES,.FALSE.)
 
      IF ( LNUL(XL) ) GOTO 99
      IF ( BONUL(XL,PAS) ) GOTO 99

      CALL SCUMW(DSREF)
      CALL SCUMR(
     >           DUM,SCUM,TCUM) 

      CALL RAZDRV(3)

      DXI=PAS
      XFE=-XE
      XFS=XS-XLIM
      IF(KP .EQ. 1)  THEN
C------- Optical elmnt undergoes new positionning. Its frame becomes the new 
C        reference frame. No exit CHANGREF needed'
        XCS = 0.D0
        YCS = 0.D0
        ALS = 0.D0

        QSHROE(1) = 'XS'
        VSHROE(1) = XCE
        QSHROE(2) = 'YS'
        VSHROE(2) = YCE
        QSHROE(3) = 'ZR'
        VSHROE(3) = ALE
        VSHROE(MSR) = 3
        QSHROS(1) = 'XS'
        VSHROS(1) = XCS
        QSHROS(2) = 'YS'
        VSHROS(2) = YCS
        QSHROS(3) = 'ZR'
        VSHROS(3) = ALS
        VSHROS(MSR) = 3

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

        IF(PAS.LE.0.D0) THEN
          TEMP=XCE
          XCE=XTEMP
          XTEMP=TEMP
          TEMP=YCE
          YCE=YTEMP
          YTEMP=TEMP
        ENDIF
        XCS=-XTEMP*CL-YTEMP*SL
        YCS=XTEMP*SL-YTEMP*CL
        ALS=-ALE

        QSHROE(1) = 'XS'
        VSHROE(1) = XCE
        QSHROE(2) = 'YS'
        VSHROE(2) = YCE
        QSHROE(3) = 'ZR'
        VSHROE(3) = ALE
        VSHROE(MSR) = 3
        QSHROS(1) = 'XS'
        VSHROS(1) = XCS
        QSHROS(2) = 'YS'
        VSHROS(2) = YCS
        QSHROS(3) = 'ZR'
        VSHROS(3) = ALS
        VSHROS(MSR) = 3

      ELSEIF(KP .EQ. 3)  THEN
C------- Optical elmnt is Z-rotated. Entrance and exit frames are 
C        tilted by (normally) half the deviation. 
C Modified, FM, Dec 05
C        YCS = -YCE
C        YCE = 0.D0
        XCS = 0.D0
        YCS = -YCE/COS(ALE)
        ALS = ALE
        CALL TRANSR(
     >              MIRROR,LDUM)
        IF(MIRROR) THEN 
          ALS = -(PI - ALE)
        ENDIF

        QSHROE(1) = 'XS'
        VSHROE(1) = XCE
        QSHROE(2) = 'YS'
        VSHROE(2) = YCE
        QSHROE(3) = 'ZR'
        VSHROE(3) = ALE
        VSHROE(MSR) = 3
        QSHROS(1) = 'XS'
        VSHROS(1) = XCS
        QSHROS(2) = 'YS'
        VSHROS(2) = YCS
        QSHROS(3) = 'ZR'
        VSHROS(3) = ALS
        VSHROS(MSR) = 3

      ELSEIF(KP .EQ. 4)  THEN
C------- Implemented for AGSMM. 
C        ALE is a delta wrt. DEV/2 
        TTA = ALE - DTTA
        ALS = TTA - DTTA
        DTTA2 = DTTA / 2.D0
        XCS = -XL * SIN(DTTA2) * SIN(DTTA2 + TTA)
        YCS = -YCE/COS(ALE) - XL * SIN(DTTA) * COS(DTTA2 + TTA)
 
C ZCE
        ZCSB = -VSHROE(4)
C PHE
        IF(VSHROE(6).NE.0.D0) THEN
          XCSB = VSHROE(5) 
          ZCSB = -VSHROE(6)
        ENDIF

C X-ROT (YAW)
C        QSHROS(1) = 'YR'
c        VSHROS(1) = Vshroe(7) 
        QSHROS(1) = 'ZS'
        VSHROS(1) = -VSHROE(6)
        QSHROS(2) = 'XS'
        VSHROS(2) = -VSHROE(5)
C Z-SHIFT
        QSHROS(3) = 'ZS'
        VSHROS(3) = -VSHROE(4)
C Z-ROT (PITCH)
        QSHROS(4) = 'YS'
        VSHROS(4) = YCS
        QSHROS(5) = 'ZR'
        VSHROS(5) = ALS
C        VSHROS(MSR) = 6
        VSHROS(MSR) = 5

      ELSEIF(KP .EQ. 5)  THEN
C X-, Y-, Z-translation, followed by 
C X-, Y-, Z-rotation wrt center of optical element at X=XL/2 

        KSR = 1
        DO WHILE (KSR .LT. NSR)



          KSR = KSR+1
        ENDDO


        XLM = XS-XE
        XTEMP=XCE
        YTEMP=YCE

        IF(PAS.LE.0.D0) CALL ENDJOB('Pgm quasex. Negative integration '
     >  //'step not supported with KPOS = ',5)

        XCS=-XTEMP
        YCS=-YTEMP

        QSHROE(1) = 'XS'
        VSHROE(1) = XCE
        QSHROE(2) = 'YS'
        VSHROE(2) = YCE
        QSHROE(3) = 'ZS'
        VSHROE(3) = ALE
        QSHROE(3) = 'ZS'
        VSHROE(3) = ZCE
        VSHROE(MSR) = 4
        QSHROS(1) = 'XS'
        VSHROS(1) = XCS
        QSHROS(2) = 'YS'
        VSHROS(2) = YCS
        QSHROS(3) = 'ZR'
        VSHROS(3) = ALS
        VSHROS(MSR) = 3

      ENDIF

      IF(PAS.LE.0.D0) THEN
        TEMP=XI
        XI=XF
        XF=TEMP
        XFE=XLIM-XS
        XFS=XE
      ENDIF
      X=XI
      XLIM=XF
 
      CALL TRANSF(QSHROE,VSHROE,QSHROS,VSHROS)

      IF(NRES .GT. 0) THEN
        WRITE(NRES,100) KP,
     >  TXTF(KP)(DEBSTR(TXTF(KP)):FINSTR(TXTF(KP))),XCE,YCE,ALE
        WRITE(NRES,FMT='(/,'' Cumulative length of optical axis = '',
     >  1P,G17.9,
     >  '' m ;  Time  (for ref. rigidity & particle) = '', 
     >  1P,G14.6,'' s '')')  SCUM*UNIT(5), TCUM
      ENDIF

      IF(KSYN.EQ.1) THEN
        CALL SRLOS1(
     >              OKPRSR,LNSR)
        IF(OKPRSR) CALL PRSR(LNSR) 
      ENDIF

 99   CONTINUE
C----- Unset wedge correction, in case would have been set by MULTIPOL, BEND, etc.
      CALL INTEG3
C----- Unset coded step
      CALL CHXC1W(I0,I0)
      CALL DEPLAW(.FALSE.,I0)

      RETURN

      ENTRY QUASE2(DTTAI,ZCEI,PHII)
      DTTA = DTTAI
      ZCE = ZCEI
      PHI = PHII
      RETURN

C100   FORMAT(/,5X,'ELEMENT  DECENTRE  PAR  RAPPORT  A  L''AXE  OPTIQUE'
C     1,/,10X,'CENTRE  DE  LA  FACE  D''ENTREE  EN  X =',1P,G10.2
C     2,' CM   Y =',G10.2,' CM   INCLINAISON =',G12.4,' RAD',/)
 100  FORMAT(/,5X,'KPOS =  ',I0,'.  ',A
     >,/,10X,'X =',1P,G12.4
     >,' CM   Y =',G12.4,' cm,  tilt  angle =',G14.6,' RAD',/)
      END
