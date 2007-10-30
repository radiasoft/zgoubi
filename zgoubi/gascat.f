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
      SUBROUTINE GASCAT(PAS,D1,IT, 
     >                            D2,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      COMMON/GASC/ AI, DEN, KGA
      COMMON/PTICUL/ AM,Q,G,TO

       PPI1=D1*0.3D00
C       AM=139.58D00
       AMB=0.511D00
       AK1=0.307075D00
       AX=14.28D00
       ZX=7.14D00
       AKX=AK1*ZX/AX
C       AI=84.25D-06
C       DEN=1.2385D-03
       EX=PPI1*PPI1+AM*AM
       EPI1=DSQRT(EX)
       TPI1=EPI1-AM
       GAMBET=PPI1/AM
       BET=PPI1/EPI1
       BET2=BET*BET
       TMAX=2.0D00*AMB*GAMBET*GAMBET
       DEL=AKX/BET2*(0.5D00*DLOG(TMAX*TMAX/AI/AI)-BET2)
       DELE=DEL*DEN*PAS
       TPI2=TPI1-DELE
       IF(TPI2 .LE. 0.D0 ) CALL KSTOP(6,IT,IEX(IT),*99)
       EPI2=TPI2+AM
       PPI2=SQRT(EPI2*EPI2-AM*AM)
       D2=PPI2/0.3D00

      RETURN
 99   RETURN 1
      END
