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
      SUBROUTINE DISPU(L, 
     >                   VAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU

      PARAMETER (MXV=60, I27 = 27, I24 = I27-3) 
      DIMENSION TAB(MXV,I27), XLOC(MXV,I24)
      SAVE IT1, IT2, NPU, XLOC

      VAL = 0.D0
      II = 0
      DO 10 JT = IT1, IT2-1, IIT
        DO 10 IT = JT+1, IT2, IIT
          DVAL2 = F(L,JT) - F(L,IT)
          VAL = VAL + DVAL2*DVAL2
          II = II + 1
 10   CONTINUE
      IF(II.GT.0) THEN 
        VAL = SQRT(VAL)/DBLE(II)
      ELSE
        VAL = 1.D10
      ENDIF
      RETURN

      ENTRY DISPU2(TAB,JJ)
      IT1 = NINT(TAB(JJ,2))
      IT2 = NINT(TAB(JJ,3))
      NPU = NINT(TAB(JJ,1)-2)
      IF(NPU.GT.I24) CALL ENDJOB
     >('Pgm disloc. Too many PU locations. Max is ',I24)
      IF(NPU.NE.999) THEN 
        DO I = 1, NPU
          XLOC(JJ,I) = TAB(JJ,2+I)        
        ENDDO
      ELSE
C All PUs upstream of IR are accounted for.
      ENDIF
      RETURN

      END
