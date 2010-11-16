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
C  Brookhaven National Laboratory               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE RDPOIS(LUN,IX,JY,X,Y,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),Y(*),B(*)
 
      DIMENSION PX(500),PB(500)
 
      READ(LUN,*,ERR=99) YMI, YMA
      DY = 2.D0*YMA/(JY-1.D0)
      Y(1)=-YMA
      DO 14 J=2,JY
 14     Y(J)= Y(J-1)+DY
 
      I=0
 11   CONTINUE
      I=I+1
      READ(LUN,101,ERR=99,END=10) IB,IB,IB,BID,PX(I), BID,BID,PB(I)
 101  FORMAT(I1,I3,I4,G15.6,2F11.5,2F12.3)
      GOTO 11
 
 10   CONTINUE
      IMA=I-1
      X(1)=PX(1)
      DX=(PX(IMA)-PX(1))/(IX-1.D0)
      DO 12 I=2,IX
 12     X(I)=X(I-1)+DX
 
      B(1)=PB(1)
      B(IX)=PB(IMA)
      I1=1
      DO 13 I=2,IX-1
        II=INTERM(X(I),PX,I1,IMA)
        B(I)=PB(II) + (PB(II+1)-PB(II))/(PX(II+1)-PX(II))
     >  *(X(I)-PX(II))
 13   CONTINUE
 
      RETURN
 99   CALL ENDJOB('ERROR DURING READ ',-99)
      END
