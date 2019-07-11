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
C  -------
      SUBROUTINE PEU(E,U,
     >                   EU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(5,3), U(6,3), EU(5,3)
      I = 1
      EU(I,1) = E(I,2) * U(I,3) - E(I,3) * U(I,2)
      EU(I,2) = E(I,3) * U(I,1) - E(I,1) * U(I,3)
      EU(I,3) = E(I,1) * U(I,2) - E(I,2) * U(I,1)
      I = 2
      EU(I,1) = E(I  ,2) * U(I-1,3) - E(I  ,3) * U(I-1,2)
     >        + E(I-1,2) * U(I  ,3) - E(I-1,3) * U(I  ,2)
      EU(I,2) = E(I  ,3) * U(I-1,1) - E(I  ,1) * U(I-1,3)
     >        + E(I-1,3) * U(I  ,1) - E(I-1,1) * U(I  ,3)
      EU(I,3) = E(I  ,1) * U(I-1,2) - E(I  ,2) * U(I-1,1)
     >        + E(I-1,1) * U(I  ,2) - E(I-1,2) * U(I  ,1)
      I = 3
      EU(I,1) = E(I  ,2) * U(I-2,3) - E(I  ,3) * U(I-2,2)
     >  + 2.D0*(E(I-1,2) * U(I-1,3) - E(I-1,3) * U(I-1,2))
     >        + E(I-2,2) * U(I  ,3) - E(I-2,3) * U(I  ,2)
      EU(I,2) = E(I  ,3) * U(I-2,1) - E(I  ,1) * U(I-2,3)
     >  + 2.D0*(E(I-1,3) * U(I-1,1) - E(I-1,1) * U(I-1,3))
     >        + E(I-2,3) * U(I  ,1) - E(I-2,1) * U(I  ,3)
      EU(I,3) = E(I  ,1) * U(I-2,2) - E(I  ,2) * U(I-2,1)
     >  + 2.D0*(E(I-1,1) * U(I-1,2) - E(I-1,2) * U(I-1,1))
     >        + E(I-2,1) * U(I  ,2) - E(I-2,2) * U(I  ,1)
      I = 4
      EU(I,1) = E(I  ,2) * U(I-3,3) - E(I  ,3) * U(I-3,2)
     >  + 3.D0*(E(I-1,2) * U(I-2,3) - E(I-1,3) * U(I-2,2))
     >  + 3.D0*(E(I-2,2) * U(I-1,3) - E(I-2,3) * U(I-1,2))
     >        + E(I-3,2) * U(I  ,3) - E(I-3,3) * U(I  ,2)
      EU(I,2) = E(I  ,3) * U(I-3,1) - E(I  ,1) * U(I-3,3)
     >  + 3.D0*(E(I-1,3) * U(I-2,1) - E(I-1,1) * U(I-2,3))
     >  + 3.D0*(E(I-2,3) * U(I-1,1) - E(I-2,1) * U(I-1,3))
     >        + E(I-3,3) * U(I  ,1) - E(I-3,1) * U(I  ,3)
      EU(I,3) = E(I  ,1) * U(I-3,2) - E(I  ,2) * U(I-3,1)
     >  + 3.D0*(E(I-1,1) * U(I-2,2) - E(I-1,2) * U(I-2,1))
     >  + 3.D0*(E(I-2,1) * U(I-1,2) - E(I-2,2) * U(I-1,1))
     >        + E(I-3,1) * U(I  ,2) - E(I-3,2) * U(I  ,1)
      RETURN
      I = 5
      EU(I,1) = E(I  ,2) * U(I-4,3) - E(I  ,3) * U(I-4,2)
     >  + 4.D0*(E(I-1,2) * U(I-3,3) - E(I-1,3) * U(I-3,2))
     >  + 6.D0*(E(I-2,2) * U(I-2,3) - E(I-2,3) * U(I-2,2))
     >  + 4.D0*(E(I-3,2) * U(I-1,3) - E(I-3,3) * U(I-1,2))
     >        + E(I-4,2) * U(I  ,3) - E(I-4,3) * U(I  ,2)
      EU(I,2) = E(I  ,3) * U(I-4,1) - E(I  ,1) * U(I-4,3)
     >  + 4.D0*(E(I-1,3) * U(I-3,1) - E(I-1,1) * U(I-3,3))
     >  + 6.D0*(E(I-2,3) * U(I-2,1) - E(I-2,1) * U(I-2,3))
     >  + 4.D0*(E(I-3,3) * U(I-1,1) - E(I-3,1) * U(I-1,3))
     >        + E(I-4,3) * U(I  ,1) - E(I-4,1) * U(I  ,3)
      EU(I,3) = E(I  ,1) * U(I-4,2) - E(I  ,2) * U(I-4,1)
     >  + 4.D0*(E(I-1,1) * U(I-3,2) - E(I-1,2) * U(I-3,1))
     >  + 6.D0*(E(I-2,1) * U(I-2,2) - E(I-2,2) * U(I-2,1))
     >  + 4.D0*(E(I-3,1) * U(I-1,2) - E(I-3,2) * U(I-1,1))
     >        + E(I-4,1) * U(I  ,2) - E(I-4,2) * U(I  ,1)

      RETURN
      END
