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
      SUBROUTINE ET2RES(lunw)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical idluni, ok, gttext, empty
      integer debstr, finstr
      character(200) txt200

      save alf1, alf2, bet1, bet2, phy, phz

      ok = idluni(
     >            nl)
      if(ok) then
        open(unit=nl,file='ETparam.res',err=91)
      else
        goto 91
      endif

      write(nres,fmt='(a)') 

      M1 = -1
      ok = GTTEXT(M1,nl,'###########',
     >                              txt200)

      ok = GTTEXT(M1,nl,'FRACTIONAL PART',
     >                              txt200)
      read(txt200(70: 89),*) phy
      read(txt200(90:106),*) phz

      ok = GTTEXT(M1,nl,'EDWARDS-TENG',
     >                              txt200)      
      read(nl,fmt='(a)',err=98,end=98) txt200
      read(txt200(70: 89),*) alf1
      read(txt200(90:106),*) alf2
      read(nl,fmt='(a)',err=98,end=98) txt200
      read(txt200(70: 89),*) bet1
      read(txt200(90:106),*) bet2

c              write(*,*) ' et2res '
c              write(*,*) ' alf 1, 2,  bet 1, 2 : '
c              write(*,*)  alf1, alf2, bet1, bet2
c                read(*,*)

      RETURN

 91   continue
      write(nres,*)      
     >'      ' //
     >' SBR MATIMC , ' // 
     >'could not open file "ETparam.res", proceeding without...'
      RETURN

 98   continue
      write(nres,*)      
      write(nres,*)      
     >'      ' //
     >' Finished copying from ETparam.res to zgoubi.res ' //
     >'upon eor or eof. '
      RETURN

      entry et2re1(
     >             F011,f012,f033,f034,phyo,phzo)
      f011 = bet1
      f012 = -alf1
      f033 = bet2
      f034 = -alf2
      phyo = phy
      phzo = phz
      return

      END
