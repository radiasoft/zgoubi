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
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE SPNMAT(ID,
     >                     SMAT,TRM, SROT,TR1,TR2,TR3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SMAT(3,3)
      INCLUDE "MAXTRA.H"
      INCLUDE "C.SPIN.H"     ! COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)

      DIMENSION DLT(3,3), PROD(3,3)

        II = 1
        IT = 1 + (ID-1)*3
        SMAT(1,II) = SF(1,IT) 
        SMAT(2,II) = SF(2,IT) 
        SMAT(3,II) = SF(3,IT) 
        II = II+1
        IT = IT+1
        SMAT(1,II) = SF(1,IT) 
        SMAT(2,II) = SF(2,IT) 
        SMAT(3,II) = SF(3,IT) 
        II = II+1
        IT = IT+1
        SMAT(1,II) = SF(1,IT) 
        SMAT(2,II) = SF(2,IT) 
        SMAT(3,II) = SF(3,IT) 

        TRM = SMAT(1,1) + SMAT(2,2) + SMAT(3,3) 
        SROT = ACOS((TRM-1.D0)/2.D0)

        CALL RAZ(DLT,3*3)
        DLT(2,3) = -1.D0 
        DLT(3,2) = +1.D0 
        PROD = MATMUL(DLT,SMAT)
        TR1 = (PROD(1,1) + PROD(2,3) + PROD(3,3))/(2.D0*SROT)
        DLT(2,3) = 0.D0 
        DLT(3,2) = 0.D0 
        DLT(1,3) = +1.D0 
        DLT(3,1) = -1.D0 
        PROD = MATMUL(DLT,SMAT)
        TR2 = (PROD(1,1) + PROD(2,3) + PROD(3,3))/(2.D0*SROT)
        DLT(1,3) = 0.D0 
        DLT(3,1) = 0.D0 
        DLT(1,2) = -1.D0 
        DLT(2,1) = +1.D0 
        PROD = MATMUL(DLT,SMAT)
        TR2 = (PROD(1,1) + PROD(2,3) + PROD(3,3))/(2.D0*SROT)
        REN = SQRT(TR1*TR1+TR2*TR2+TR3*TR3)
        TR1 = TR1 / REN ; TR2 = TR2 / REN ; TR3 = TR3 / REN 

      RETURN
      END
