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
C  Brookhaven National Laboratory               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE SOLBAX(XL,BOSQ,RO2,X,
     >                                BX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BX(*)
                
      X2 = X*X
      U = X2 + RO2
      U2 = U*U
      U3 = U2*U
      SQU = SQRT(U)
      DU = 2.D0 * X
      XL2 = (X-XL)*(X-XL)
      V = XL2 + RO2
      V2 = V*V
      V3 = V2*V
      SQV = SQRT(V)
      DV = 2.D0 * (X-XL)
C------- BX(n)=-d(n-1)B/dX^(n-1)
C        BX(1)=BOSQ * ( X/SQU - (X-XL)/SQV )
        BX(1)=BOSQ * (
     >x/Sqrt(ro2 + x**2) + (-x + xl)/Sqrt(ro2 + (x - xl)**2)
     > )

        BX(2)=BOSQ * (
     >ro2*((ro2 + x**2)**(-1.5) - (ro2 + (x - xl)**2)**(-1.5))
     > )

        BX(3)=BOSQ * (
     >ro2*((-3*x)/(ro2 + x**2)**2.5 + 
     >(3*(x - xl))/(ro2 + (x - xl)**2)**2.5)
     > )

        BX(4)=BOSQ * (
     >  3*ro2*(4/(ro2 + x**2)**2.5 + 
     >5*ro2*(-(ro2 + x**2)**(-3.5) + (ro2 + (x - xl)**2)**(-3.5)) - 
     >4/(ro2 + (x - xl)**2)**2.5)
     > )

        BX(5)=BOSQ * (
     >   3*ro2*((-20*x)/(ro2 + x**2)**3.5 + 
     >    (20*(x - xl))/(ro2 + (x - xl)**2)**3.5 + 
     >    5*ro2*((7*x)/(ro2 + x**2)**4.5 + 
     >       (-7*x + 7*xl)/(ro2 + (x - xl)**2)**4.5))
     > )

        BX(6)= BOSQ * (
     >45*ro2*(8/(ro2 + x**2)**3.5 + 
     >21*ro2**2*((ro2 + x**2)**(-5.5) - (ro2 + (x - xl)**2)**(-5.5)) + 
     >28*ro2*(-(ro2 + x**2)**(-4.5) + (ro2 + (x - xl)**2)**(-4.5)) - 
     >8/(ro2 + (x - xl)**2)**3.5)
     > )

      RETURN
      END
