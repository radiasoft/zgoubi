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
C  Brookhaven National Laboratory                                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE CALTRI(N,S,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C----------------------------------------------------------------------------
C
C     N   : MODE   1 RECHERCHE ECHELLE MAX
C           MODE   2 TRACE
C     S   : POSITION X
C     Y   : POSITION Y 
C
C----------------------------------------------------------------------------
      COMMON/VEN/Y1,S1,D,DQ,DS,DM,DSS,DPU,TETA1,
     >           SMIN,SMAX,YMIN,YMAX,SB0,YB0
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

C     SAVE POSITION IN COMMON FOR FUTUR USE

      S1=S
      Y1=Y
C
C     ADD OFFSET
C
      S2=S1+SB0
      Y2=Y1+YB0
      IF(N.LT.2) THEN
C
C        MAX SEARCH (NO PLOT)
C
         IF(S2.GT.SMAX) SMAX=S2
         IF(S2.LT.SMIN) SMIN=S2
         IF(Y2.GT.YMAX) YMAX=Y2
         IF(OKECH) YMAX=YMA
         IF(Y2.LT.YMIN) YMIN=Y2
         IF(OKECH) YMIN=YMI
      ELSE
C
C        PLOT
C
         CALL VECTPL(S2,Y2,2)

      ENDIF
      RETURN
      END  
