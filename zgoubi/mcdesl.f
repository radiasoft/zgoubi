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
      SUBROUTINE MCDESL(DL,I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     *******************************
C     APPELE PAR ESL , SEPARATEUR...
C     *******************************
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"      ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE "MAXTRA.H"
      INCLUDE "C.DESIN.H"     ! COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
C     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
C     >,AMS ,AMP,ENSTAR,BSTAR,TDVM,TETPHI(2,MXT)
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AMASS,Q,G,TO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)

C      DO 3 IT=1,IMAX
      DO 3 IT=I1,I2

C------- IEX<-1<=> PARTICULE STOPPEE
        IF(IEX(IT) .LT. -1) GOTO 3

        Y = F(2,IT)
        T = F(3,IT) * UNIT(2)
        Z = F(4,IT)
        P = F(5,IT) * UNIT(4)
        SAR= F(6,IT)
        TAR= F(7,IT)
        DP= F(1,IT)
        QBR=Q*BORO*DP
        AMT = AMQ(1,IT)
        QT = AMQ(2,IT)
        CALL MCDES(DL,IEX(IT),Y,T,Z,P,DL,QBR,SAR,TAR,IT,AMT,Q,BORO,UN,1)
        F(2,IT) = Y
        F(3,IT) = T / UNIT(2)
        F(4,IT) = Z
        F(5,IT) = P / UNIT(4)
        F(6,IT) = SAR
        F(7,IT) = TAR
        DP=QBR/(Q*BORO)
        F(1,IT) = DP
        AMQ(1,IT) = AMT
        AMQ(2,IT) = QT

 3    CONTINUE

      RETURN
      END
