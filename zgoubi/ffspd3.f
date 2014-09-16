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

C   T=TTA, TTARF=omega
      FUNCTION FFSPD3(XB,YB,T,TTARF,a,b,c)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      COST=COS(TTARF+T)
      SINT=SIN(TTARF+T)
      FFSPD3= 
     - (XB - (-c/(a*COST+b*SINT))*COST)  * (-b) +
     - (YB - (-c/(a*COST+b*SINT))*SINT)  * (a)



!      COST=COS(TTARF+T)
!      SINT=SIN(TTARF+T)
!      FFSPD3= 
!     - (XB - (-c/(a*COST+b*SINT))*COST)  *
!     - ((c/(a*COST+b*SINT)**2 *(-a*SINT+b*COST)*COST)  + 
!     - (c/(a*COST+b*SINT))*SINT) +
!     - (YB - (-c/(a*COST+b*SINT))*SINT)  * 
!     - ((c/(a*COST+b*SINT)**2 *(-a*SINT+b*COST)*SINT)  + 
!     - (-c/(a*COST+b*SINT))*COST)
!        WRITE(*,*) TTARF,T
      RETURN

      ENTRY FFSPDD3(XB,YB,T,TTARF,a,b,c)

       FFSPDD3= 
     - ((-c/(a*COST+b*SINT)**2 *(-a*SINT+b*COST)*COST) * (-b) +
     - (c/(a*COST+b*SINT))*(-SINT)) * (-b) +
     - (-c/(a*COST+b*SINT)**2 *(-a*SINT+b*COST)*SINT) * (a)   + 
     - (c/(a*COST+b*SINT))*(COST) * (a)  



!       FFSPDD3= 
!     - ((-c/(a*COST+b*SINT)**2 *(-a*SINT+b*COST)*COST) +
!     - (c/(a*COST+b*SINT))*(-SINT)) *
!     - ((c/(a*COST+b*SINT)**2 *(-a*SINT+b*COST)*COST)  + 
!     - (c/(a*COST+b*SINT))*SINT)  +
!     - (XB - (-c/(a*COST+b*SINT))*COST)  *
!     - ((-2*c/(a*COST+b*SINT)**3 *(-a*SINT+b*COST)**2 * COST)*
!     - (-c/(a*COST+b*SINT)**2 *(-a*SINT+b*COST)*SINT) + 
!     - (c/(a*cos(T)+b*sin(T)))*(COST))  
      RETURN
      END
