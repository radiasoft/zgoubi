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
      SUBROUTINE KRTAX(X1,BX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BX(*)
      INCLUDE 'PARIZ.H'
      INCLUDE "XYZHC.H"
C      COMMON//XH(MXX),YH(MXY),ZH(IZ),HC(ID,MXX,MXY,IZ),IXMA,JYMA,KZMA
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      COMMON/RIGID/ BORO,DPREF,DP,BR
 
C      REAL*8 C(6),F(5)
      DIMENSION C(6),F(5)
 
C... B(X) ==  HC(ID,IX,1,1),  X(IX) == XH(IX)
C... X(2)- X(1) == LONG. D'UNE MAILLE
      DX=XH(2)-XH(1)
      DX2=DX*DX
 
C... X1 == X COURANT
      IAC=(X1-XH(1))/DX+1.5D0
 
C... GRILLE 1-D A 5 POINTS OU 2-D A 5*5 POINTS
      IAC=MAX0(IAC,3)
      IAC=MIN0(IAC,IXMA-2)
 
C.... DX == X/COURANT - X/CENTREGRILLE  LE TOUT /DX
      A=(X1-XH(IAC))/DX
 
      F(1)=HC(ID,IAC-2,1,1)/BR
      F(2)=HC(ID,IAC-1,1,1)/BR
      F(3)=HC(ID,IAC  ,1,1)/BR
      F(4)=HC(ID,IAC+1,1,1)/BR
      F(5)=HC(ID,IAC+2,1,1)/BR
 
      E1 =     F(1)+F(2)+F(3)+F(4)+    F(5)
      E2 = -2.D0*F(1)-F(2)     +F(4)+ 2.D0*F(5)
      E3 =  4.D0*F(1)+F(2)     +F(4)+ 4.D0*F(5)
      E4 = -8.D0*F(1)-F(2)     +F(4)+ 8.D0*F(5)
      E5 = 16.D0*F(1)+F(2)     +F(4)+16.D0*F(5)
 
      C(5)=(35.D0*E5-155.D0*E3+ 72.D0*E1         )/288D0
      C(3)=(            E3-  2.D0*E1-62.D0*C(5))/14d0
      C(1)=(                    E1-34.D0*C(5)-10.D0*C(3))/5D0
 
      C(4)=(10.D0*E4-34.D0*E2         )/144D0
      C(2)=(           E2-34.D0*C(4))/10D0
 
      C(6)=0D0
 
C      CALL DRVG(4,C,A,G,DG,D2G,D3G,D4G,D5G,D6G)
C      A=A*DX
C      DF  =(F(4)-F(2))/2/DX
C      D2F =(-F(5)+16*F(4)-30*F(3)+16*F(2)-F(1))/12/DX2
C      D3F =(F(5)-2*F(4)+2*F(2)-F(1))/2/DX2/DX
C      D4F =(F(5)-4*F(4)+6*F(3)-4*F(2)+F(1))/DX2/DX2
C      C(1)=F(3)
C      C(2)=DF
C      C(3)=D2F/2
C      C(4)=D3F/6
C      C(5)=D4F/24
 
C----- BX(n+1) = dnBX/dXn
      BX(1)=  C(1)+(C(2)+(C(3)+(C(4)+(C(5)+C(6)*A)*A)*A)*A)*A
      BX(2)= (C(2)+(2*C(3)+(3*C(4)+(4*C(5)+5*C(6)*A)*A)*A)*A)/DX
      BX(3)= (2.D0*C(3)+(6.D0*C(4) +(12.D0*C(5) +20.D0*C(6)*A)*A)*A)/DX2
      BX(4)= (6.D0*C(4)+(24.D0*C(5) + 60.D0*C(6)*A)*A)/DX/DX2
      BX(5)= (24.D0*C(5)+120.D0*C(6)*A)/DX2/DX2
      BX(6)=0D0
 
      RETURN
      END
