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

CCCCCC using newton's method this procedure will solve the equation: r=RM*exp(b(r)*theta) where b(r) is defined earlier and theta=fixed here but varying in rtnewt2.

      SUBROUTINE RSOLVE(RTNEWT2,D,sol0,RS,FA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
C     Initial guess
      RS=sol0
      FRLIM=0.001
      RACC=0.001D0
      JMAX=10000
      thet=RTNEWT2
      
      DO 11 J=1,JMAX
         fr = funct(RS,thet,D,FA)
         frprime = functd(RS,thet,D,FA)
C        write(*,*) DX
C     DEB RAJOUT
        
         IF (ABS(fr).GT.FRLIM) THEN
            deltar = fr/frprime
            RS = RS - deltar            
         ENDIF


         IF(ABS(fr).LT.RACC) RETURN 
 11      CONTINUE
C      PAUSE 'rtnewt2 exceeding maximum iterations'

      RETURN
      END
