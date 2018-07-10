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
      SUBROUTINE SOLBAX(XL,BO2,RO2,X,
     >                               BX)
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
C----- BX(n)=-d(n-1)B/dX^(n-1)
C      BX(1)=BO/2 * ( X/SQU - (X-XL)/SQV )
      BX(1)=BO2 * (
     >X/SQRT(RO2 + X**2) + (-X + XL)/SQRT(RO2 + (X - XL)**2) )

      BX(2)=BO2 * (
     >RO2*((RO2 + X**2)**(-1.5) - (RO2 + (X - XL)**2)**(-1.5)) )

      BX(3)=BO2 * (
     >RO2*((-3.D0*X)/(RO2 + X**2)**2.5 + 
     >(3.D0*(X - XL))/(RO2 + (X - XL)**2)**2.5) )

      BX(4)=BO2 * (
     >3.D0*RO2*(4.D0/(RO2 + X**2)**2.5 + 
     >5.D0*RO2*(-(RO2 + X**2)**(-3.5) + (RO2 + (X - XL)**2)**(-3.5)) - 
     >4.D0/(RO2 + (X - XL)**2)**2.5) )

C              write(*,*) ' SOLBA  BO2, BX(4) = ', BO2,BX(4)
c              write(*,*) ' SOLBA  BX(5) = ',BX(5),BX(6),BX(7),BX(8)
C     >,BX(9),BX(10),BX(11),BX(12),BX(13),BX(14),BX(15)


      BX(5)=BO2 * (
     > 3.D0*RO2*((-20.D0*X)/(RO2 + X**2)**3.5 + 
     >  (20.D0*(X - XL))/(RO2 + (X - XL)**2)**3.5 + 
     >  5.D0*RO2*((7.D0*X)/(RO2 + X**2)**4.5 + 
     >     (-7*X + 7.D0*XL)/(RO2 + (X - XL)**2)**4.5)) )
 
      BX(6)= BO2 * (
     >45.D0*RO2*(8.D0/(RO2 + X**2)**3.5 + 
     >21.D0*RO2**2*((RO2 + X**2)**(-5.5) -(RO2 + (X - XL)**2)**(-5.5))+ 
     >28.D0*RO2*(-(RO2 + X**2)**(-4.5) + (RO2 + (X - XL)**2)**(-4.5)) - 
     >8.D0/(RO2 + (X - XL)**2)**3.5) )

      RETURN
      END
