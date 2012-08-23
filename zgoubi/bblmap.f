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
        subroutine bblmap(linearmap,tunex,tuney,tunez,
     >              alphax,betax,alphay,betay,alfmom,ga,bet,clight,circ)
C S. White & F. Meot, Jan. 2012
        implicit double precision (a-h,o-z)
        double precision linearmap(6,6) 

        twopi = 4.d0*dasin(1.0d0)
        gx = (1.d0+alphax*alphax)/betax
	cx=dcos(twopi*tunex)
	sx=dsin(twopi*tunex)
        linearmap(1,1) = cx + alphax*sx
        linearmap(2,1) = -gx*sx
        linearmap(1,2) = betax*sx
        linearmap(2,2) = cx - alphax*sx
        gy = (1+alphay*alphay)/betay
	cy=dcos(twopi*tuney)
	sy=dsin(twopi*tuney)
        linearmap(3,3) = cy + alphay*sy
        linearmap(4,3) = -gy*sy
        linearmap(3,4) = betay*sy
        linearmap(4,4) = cy - alphay*sy
	cz=dcos(twopi*tunez)
	sz=dsin(twopi*tunez)
        szspz = (alfmom-1/(ga*ga))*bet*clight/
     >            (tunez*bet*clight*twopi/circ)
        linearmap(5,5) = cz
        linearmap(6,5) = sz/szspz
        linearmap(5,6) = -sz*szspz
        linearmap(6,6) = cz

!         Set coupling terms to 0

        linearmap(1,3)=0.d0
        linearmap(1,4)=0.d0
        linearmap(1,5)=0.d0
        linearmap(1,6)=0.d0
        linearmap(2,3)=0.d0
        linearmap(2,4)=0.d0
        linearmap(2,5)=0.d0
        linearmap(2,6)=0.d0

        linearmap(3,1)=0.d0
        linearmap(3,2)=0.d0
        linearmap(3,5)=0.d0
        linearmap(3,6)=0.d0
        linearmap(4,1)=0.d0
        linearmap(4,2)=0.d0
        linearmap(4,5)=0.d0
        linearmap(4,6)=0.d0

        linearmap(5,1)=0.d0
        linearmap(5,2)=0.d0
        linearmap(5,3)=0.d0
        linearmap(5,4)=0.d0
        linearmap(6,1)=0.d0
        linearmap(6,2)=0.d0
        linearmap(6,3)=0.d0
        linearmap(6,4)=0.d0

        return
        end
