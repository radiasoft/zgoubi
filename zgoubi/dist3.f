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
      SUBROUTINE DIST3(L,NOELA,NOELB,
     >                               VAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      
      DIMENSION XYZ1(MXJ), XYZ2(MXJ)
C Get 1st and last PU #s, as well as corresponding NOELs, for PUs located in range noela-noelb
      CALL PCKUP7(NOELA,NOELB, 
     >                         IPUI,IPUF,NOELI,NOELF)
      VAL = 0.D0
      KPU = 0
      DO IPU = IPUI, IPUF-1
        CALL PCKUP5(IPU
     >                 ,XYZ1)
        DO JPU = IPU+1, IPUF
          KPU = KPU+1
          CALL PCKUP5(JPU
     >                   ,XYZ2)
          VAL = VAL + ABS(XYZ1(L)-XYZ2(L))
        ENDDO
      ENDDO

      RETURN
      END
