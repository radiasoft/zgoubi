C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <meot@lpsc.in2p3.fr>
C  Service Acc�lerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE ESLVIR(IOPT,YL,ZL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------
C     SECTION SANS Champ
C     ------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
     > ,YCH,ZCH
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
     1,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/RIGID/ BORO,DPREF,DP,BR
 
      DLY=YL
      DLZ=ZL
      IF(IOPT .EQ. 1) THEN
        DLY = A(NOEL,1)
        DLZ = A(NOEL,2)
        IF(NRES.GT.0) WRITE(NRES,109) DLY,DLZ
 109    FORMAT(/,10X,'ESPACE  VIRTUEL Y,Z=',
     S         F12.5,1X,F12.5,'  CM',/)
      ENDIF
 
 
      DO 1 I=1,IMAX
C-------- IEX <-1 <=> PARTICULE STOPPEE
         IF(IEX(I) .LT. -1) GOTO 1
 
         FO(2,I)=FO(2,I)+DLY*TAN(FO(3,I)*1.D-3)
         F(2,I)=FO(2,I)
         FO(4,I)=FO(4,I)+DLZ*TAN(FO(5,I)*1.D-3)
         F(4,I)=FO(4,I)
 1    CONTINUE
 
 
      IF(NRES.GT.0) WRITE(NRES,101) IEX(1),(F(J,1),J=1,7)
  101 FORMAT(' TRAJ 1 IEX,D,Y,T,Z,P,S,time :',I3,1P,5G12.4,2G17.5)

      RETURN 
      END
