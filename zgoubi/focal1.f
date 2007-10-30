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
      SUBROUTINE FOCAL1(IRF,MX1,MX2,XI,YI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------------
C     CALCULE LA POSITION XI DU POINT DE FOCALISATION
C     APPELE PAR SPGM MATRIX.
C     -----------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
 
      DIMENSION TTI(MXT)
 
      TTI(1)=TAN(F(3,IRF)*.001D0)
      STY=TTI(1)*F(2,IRF)
      ST2=TTI(1)*TTI(1)
      ST=TTI(1)
      SY=F(2,IRF)
      MXI=1
      DO 1 I=MX1,MX2
         IF(I .EQ. IRF) GOTO 1
         MXI=MXI+1
         TTI(I)=TAN(F(3,I)*1.D-3)
         STY=STY+TTI(I)*F(2,I)
         ST2=ST2+TTI(I)*TTI(I)
         ST=ST+TTI(I)
         SY=SY+F(2,I)
 1    CONTINUE
      XI=-(STY-(ST*SY)/MXI)/(ST2-(ST*ST)/MXI)
      YI =F(2,IRF)+XI*TTI(1)
 
      RETURN
      END
