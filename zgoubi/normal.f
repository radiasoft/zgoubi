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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE NORMAL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ****************************************************************
C     NORMALLY FOR MULTITURN TRACKING:
C     - RENORMALISES PHASE-SPACE COORDINATES ( => EMITTANCE=CONSTANT )
C     - RENORMALISES SPIN COORDINATES ( => SPIN=1 )
C     ****************************************************************
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
C      COMMON/SYMPL/ NRME,NRMS
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
 
      DIMENSION EMY(MXT),EMZ(MXT)
      SAVE EMY, EMZ
 
C     ... TRANSVERSE PHASE_SPACES
      NRME =A(NOEL,1)
C     ... BETA, ALPHA
      BY   =A(NOEL,10)
      AY   =A(NOEL,11)
C     ... BETATRON CLOSED ORBIT
      YO   =A(NOEL,12)
      YPO  =A(NOEL,13)
C     ...RG, RG'
      DY   =A(NOEL,14)
      DPY  =A(NOEL,15)
 
      BZ   =A(NOEL,20)
      AZ   =A(NOEL,21)
      ZO   =A(NOEL,22)
      ZPO  =A(NOEL,23)
 
C     ... SPIN
      NRMS=A(NOEL,30)
 
      IF(NRME .NE. 1) NRME=0
      IF(NRMS .NE. 1) NRMS=0
      IF(NRME+NRMS .EQ. 0) WRITE(NRES,100)
 100  FORMAT(/,10X,'NO  RENORMALISATION  REQUIRED')
 
      IF(NRME .EQ. 1) THEN
 
        GY=(1.D0+AY*AY)/BY
        GZ=(1.D0+AZ*AZ)/BZ
        IF(NRES.GT.0) WRITE(NRES,110) IMAX,BY,BZ,AY,AZ,GY,GZ
     >  ,YO,ZO,YPO,ZPO,DY,DPY
 110    FORMAT(/,10X,'COURANT  INVARIANT OF  THE ',I3,'  PARTICLES',2X
     >  ,'RENORMALIZED  TO  INITIAL  VALUES  AT   EACH  PASS.'
     >  ,//,10X,'TWISS  PARAMETERS:'
     >  ,/,20X,'BETA-Y  =',1P,G12.4,' M         BETA-Z  =',G12.4,' M'
     >  ,/,20X,'ALPHA-Y =',   G12.4,'           ALPHA-Z =',G12.4
     >  ,/,20X,'GAMMA-Y =',1P,G12.4,' /M        GAMMA-Z =',G12.4,' /M'
     >  ,/,20X,'Y O.F.  =',   G12.4,' M         Z O.F.  =',G12.4,' M'
     >  ,/,20X,'Y''O.F.  =',  G12.4,' RD        Z''O.F. =',G12.4,' RD'
     >  ,/,20X,'RG      =',   G12.4,' M '
     >  ,/,20X,'RG''     =',  G12.4,//)
 
        IF(IPASS.GT.1) THEN
          DO 10 IT=1,IMAX
C           ... HORIZONTAL
            Y =F(2,IT)*1.D-2 -YO  -DY *(1.D0-FO(1,IT) )
            YP=F(3,IT)*1.D-3 -YPO -DPY*(1.D0-FO(1,IT) )
            AN=SQRT(EMY(IT) / (GY*Y*Y+2.D0*AY*Y*YP+BY*YP*YP) )
            F(2,IT) = ( Y*AN   +YO  +DY *(1.D0-FO(1,IT) ))*1.D2
            F(3,IT) = ( YP*AN  +YPO +DPY*(1.D0-FO(1,IT) ))*1.D3
C           ... VERTICAL
            Y =F(4,IT)*1.D-2 -ZO
            YP=F(5,IT)*1.D-3 -ZPO
            AN=SQRT(EMZ(IT) / (GZ*Y*Y+2.D0*AZ*Y*YP+BZ*YP*YP) )
            F(4,IT) = ( Y*AN   +ZO )*1.D2
            F(5,IT) = ( YP*AN  +ZPO)*1.D3
 
 10       CONTINUE
        ELSEIF(IPASS .EQ. 1) THEN
          DO 11 IT=1,IMAX
            Y =FO(2,IT)*1.D-2 -YO  -DY *(1.D0-FO(1,IT) )
            YP=FO(3,IT)*1.D-3 -YPO -DPY*(1.D0-FO(1,IT) )
            EMY(IT)=GY*Y*Y+2.D0*AY*Y*YP+BY*YP*YP
            Y =FO(4,IT)*1.D-2 -ZO
            YP=FO(5,IT)*1.D-3 -ZPO
            EMZ(IT)=GZ*Y*Y+2.D0*AZ*Y*YP+BZ*YP*YP
 11       CONTINUE
        ENDIF
 
        IF(NRES.GT.0) WRITE(NRES,111) (I,EMY(I),EMZ(I),I=1,IMAX)
 111    FORMAT(5X,'PARTICLE # ',I3,3X,'EY/PI =',1P,G12.4,' M.RD'
     >  ,5X,'EZ/PI =',1P,G12.4,' M.RD')
 
      ENDIF
 
      IF(NRMS .EQ. 1) THEN
 
        IF(NRES.GT.0) WRITE(NRES,120) IMAX
 120    FORMAT(/,10X,80('.'),//,10X,'SPIN  OF  THE ',I3,'  PARTICLES'
     >  ,2X,'RENORMALIZED  TO  INITIAL  VALUES  AT   EACH  PASS')
 
        DO 20 IT=1,IMAX
          AN=SQRT(SF(1,IT)*SF(1,IT)+SF(2,IT)*SF(2,IT)+SF(3,IT)*SF(3,IT))
          DO 20 J=1,3
            SF(J,IT)=SF(J,IT)/AN*SI(4,IT)
 20     CONTINUE
 
      ENDIF
 
      RETURN
      END
