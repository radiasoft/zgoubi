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
C  François Meot <fmeot@bnl.gov>
C  Brookhaven National Laboratory 
C  C-AD, Bldg 911
C  Upton, NY, 11973, USA
C  -------
      FUNCTION AGSMMA(DEV)
C     ------------------------------------------
C     Convert rigidity to amperes in main coils.
C     R. Thern, E. Blesser, AGS/AD/Tech.Note 424
C     ------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      PARAMETER (AM4=-1.070694D3, AM3=+5.836957D2, AM2=-1.119693D2, 
     >AM1=+1.282979D0, A0=8.224488D0, A1=-3.408321D-4, A2=+7.313631D-6, 
     >A3=-5.642907D-8,A4=+1.966953D-10,A5=-3.143464D-13,A6=1.897916D-16)

      BL = BORO*1.D-3 *(DPREF+HDPRF) * DEV
      BL1 = 1.D0 / BL
      AGSMMA = AM1 +(AM2 + (AM3 + AM4*BL1)*BL1)*BL1 + 
     >(A0 + (A1 + (A2 + (A3 + (A4 + (A5 + A6*BL)*BL)*BL)*BL)*BL)*BL)*BL
      RETURN
      END
