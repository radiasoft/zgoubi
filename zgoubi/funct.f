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
      FUNCTION FUNCT(R,A,RM,FA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.SPIRALE.H"     ! COMMON/spiral_ent/UMEG,ASP0,ASP1,ASP2,ASP3 
      INCLUDE "C.SPIRALX.H"     ! COMMON/spiral_ext/UMEGs,ASPS0,ASPS1,ASPS2,ASPS3
      SAVE XI,EBT
      IF (FA==1.0D0) THEN
        XI= ASP0+ASP1*R+ASP2*R**2+ASP3*R**3
      ELSE
        XI= ASPS0+ASPS1*R+ASPS2*R**2+ASPS3*R**3  
      ENDIF
!      WRITE(*,*) XI
      EBT=EXP(A/TAN(XI))
      FUNCT= R-RM*EBT

      RETURN

      ENTRY FUNCTD(R,A,RM,FA)
      IF (FA==1.0D0) THEN
        XI1=ASP1 +2.D0*ASP2*R+3.D0*ASP3*R**2
      ELSE
        XI1=ASPS1+2.D0*ASPS2*R+3.D0*ASPS3*R**2    
      ENDIF
!      WRITE(*,*) XI1
      FUNCTD= 1.D0+RM*A*(1.D0+1.D0/(TAN(XI))**2)*EBT*XI1

      RETURN
      END
