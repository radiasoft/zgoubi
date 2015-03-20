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
      FUNCTION FUNCT(r,a,RM,FA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.SPIRALE.H"     ! COMMON/spiral_ent/UMEG,ASP0,ASP1,ASP2,ASP3 
      INCLUDE "C.SPIRALX.H"     ! COMMON/spiral_ext/UMEGs,ASPS0,ASPS1,ASPS2,ASPS3
      SAVE xi,EBT
      IF (FA==1.0) THEN
      xi= ASP0+ASP1*r+ASP2*r**2+ASP3*r**3
      ELSE
      xi= ASPS0+ASPS1*r+ASPS2*r**2+ASPS3*r**3  
      ENDIF
!      write(*,*) xi
      EBT=exp(a/tan(xi))
      FUNCT= r-RM*EBT

      RETURN


      ENTRY FUNCTD(r,a,RM,FA)
      IF (FA==1.0) THEN
      xi1=ASP1+2*ASP2*r+3*ASP3*r**2
      ELSE
      xi1=ASPS1+2*ASPS2*r+3*ASPS3*r**2    
      ENDIF
!      write(*,*) xi1
      FUNCTD= 1+RM*a*(1+1/(tan(xi))**2)*EBT*xi1

      RETURN
      END
