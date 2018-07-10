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
C  Upton, NY, 11973,  USA
C  -------
      SUBROUTINE YMOINY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
 
      IF(NRES.GT.0) WRITE(NRES,100)
100   FORMAT(/,40X,'*****************   YMOINSY   ******************',/)
      DO 1 I=1,IMAX
         DO 1 J=2,5
            F(J,I)=-F(J,I)
1     CONTINUE

      IF(NRES.GT.0) THEN
        WRITE(NRES,101) IEX(1),-1.D0+F(1,1),(F(J,1),J=2,7)
 101    FORMAT('TRAJ #1 IEX,D,Y,T,Z,P,S,time :',
     >  I3,1P,5E14.6,1X,E15.7,1X,E13.5)
        IF(KSPN.EQ.1) WRITE(NRES,102) IEX(1),(SF(I,1),I=1,4)
 102    FORMAT('TRAJ #1 SX, SY, SZ, |S| :',1X,I2,2X,1P,4(E14.6,1X))
      ENDIF

      RETURN
      END
