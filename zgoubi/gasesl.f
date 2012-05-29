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
      SUBROUTINE GASESL(DL,I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -------------------------------
C     Called by ESL, SEPARA ...
C     -------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI

C      DO 3 IT=1,IMAX
      DO 3 IT=I1,I2
C------- IEX<-1<=> PARTICULE STOPPEE
        IF(IEX(IT) .LT. -1) GOTO 3
 
        Y = F(2,IT)
        T = F(3,IT)*.001D0
        Z = F(4,IT)
        P = F(5,IT)*.001D0
        SAR= F(6,IT)
        DP= F(1,IT)
        QBR0=Q*BORO*DP
        CALL GASCAT(DL,QBR0,IT,
     >                        QBR,*3)
        F(2,IT) = Y
        F(3,IT) = T*1.D3
        F(4,IT) = Z
        F(5,IT) = P*1.D3
        F(6,IT) = SAR
        DP=QBR/(Q*BORO)
        F(1,IT) = DP
 
 3    CONTINUE
 
      END
