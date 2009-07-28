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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE BAR2XY(BZAR,R,IDB,
     >                              BZ0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-------------------------------------------------
C Convert BZ(R,A) from Z_axis-polar into BZ(X,Y,Z)
C-------------------------------------------------
      PARAMETER (MDA=5,MDR=5)
      DIMENSION BZAR(MDA,MDR)
      DIMENSION BZ0(5,5)

C----- Field BZ and derivatives wrt R,A at current position

C------- Transformation polar-cartesian. BZ and X,Y,Z deriavtives 
C     computed from BZ and derivatives wrt A,R at current position
         R1=1.D0/R
         R2=R1*R1
         BZ   =   BZAR(1,1)
         BZX   =   BZAR(2,1)*R1
         BZX4  =(((BZAR(5,1) -8.D0*BZAR(3,1))*R1+6.D0*BZAR(3,2)-
     >3.D0*BZAR(1,2))*R1+ 3.D0*BZAR(1,3))*R2
         BZX3Y =(( 6.D0*BZAR(2,1)-3.D0*BZAR(4,1)*R1+BZAR(4,2)-
     >8.D0*BZAR(2,2))*R1+ 3.D0*BZAR(2,3))*R2
         BZX2Y2=(((6.D0*BZAR(3,1)*R1-4.D0*BZAR(3,2)+2.D0*BZAR(1,2))*R1-
     >    2.D0*BZAR(1,3)+BZAR(3,3))*R1+BZAR(1,4))*R1
         BZXY3 =((6.D0*(BZAR(2,2)-BZAR(2,1))*R1-3.D0*BZAR(2,3))*R1+
     >BZAR(2,4))*R1
         BZXXX = ( BZAR(4,1)*R1 + 3.D0*BZAR(2,2) - 2.D0*BZAR(2,1) )*R2
         BZXXY =(( BZAR(3,2) - 2.D0*BZAR(3,1)*R1 - BZAR(1,2) )*R1 + 
     >BZAR(1,3))*R1
         BZXYY =   BZAR(2,3)*R1 + 2.D0*( BZAR(2,1) - BZAR(2,2) )*R2    
         BZXX  = ( BZAR(3,1)*R1+BZAR(1,2))*R1
         BZXY  = ( BZAR(2,2)-BZAR(2,1))*R1
 
C Reminicences from ancient times
      BZ0(1,1)=BZ
      BZ0(2,1)=BZX
C      BZ0(1,2)=BZY
      BZ0(3,1)=BZXX
      BZ0(2,2)=BZXY
C      BZ0(1,3)=BZYY
      IF (IDB.GE.3) THEN
         BZ0(4,1)=BZXXX
         BZ0(3,2)=BZXXY
         BZ0(2,3)=BZXYY
C        BZ0(1,4)=BZYYY
         IF (IDB.GE.4) THEN
            BZ0(5,1)=BZX4
            BZ0(4,2)=BZX3Y
            BZ0(3,3)=BZX2Y2
            BZ0(2,4)=BZXY3
C           BZ0(1,5)=BZY4
         ENDIF
      ENDIF
      RETURN
      END
