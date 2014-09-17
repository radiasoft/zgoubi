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
C  Upton, NY, 11973
C  -------
      SUBROUTINE IMPMOD
     >(PGMNAM,NRES,OKCPLD,F011,f012,f033,f034,phy,phz,Cstrn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) PGMNAM
      LOGICAL OKCPLD, OKCPS

      SAVE OKCPS, F011S, F012S, F033S, F034S, PHYS, PHZS, CSTRNS
      logical first, idluni
      save first
         data first / .true. /

      DATA OKCPS / .FALSE. /

         if(first) then
           first = idluni(
     >                lunmod)
           first = .false.           
            open(lunmod,file='zgoubi.IMPMOD.out')
            write(lunmod,*) '# bet1, alf1, bet2, alf2, q1, q2, |C|'
            write(lunmod,*) '# '
         endif


      OKCPS=OKCPLD
      F011S=F011
      f012S=f012
      f033S=f033
      f034S=f034
      phyS=phy
      phzS=phz
      CstrnS=Cstrn

      ENTRY IMPMO1

      IF(NRES .GT. 0) then
       IF(OKCPS) THEN
        WRITE(NRES,*)
        WRITE(NRES,*) '--------------------------------------'
        WRITE(NRES,*) ' Pgm impmod, called by ',pgmnam
        WRITE(NRES,*) ' Coupled optics assumed.  '
        WRITE(NRES,*)
        WRITE(NRES,*) ' Coupled modes : '
        WRITE(NRES,*) ' bet1, alf1 :    ',    F011,-f012
        WRITE(NRES,*) ' bet2, alf2 :    ',    f033,-f034
        WRITE(NRES,*) ' Q1, Q2 :        ',    phy,phz  
        WRITE(NRES,*) ' Coupling strength :     ',    Cstrn
        WRITE(NRES,*)
        WRITE(NRES,*) '--------------------------------------'
        WRITE(NRES,*)
       ELSE
        WRITE(NRES,*)
        WRITE(NRES,*) '--------------------------------------'
        WRITE(NRES,*) ' Pgm impmod, called by ',pgmnam
        WRITE(NRES,*) ' Un-coupled optics assumed.  '
        WRITE(NRES,*)
        WRITE(NRES,*) '--------------------------------------'
        WRITE(NRES,*)
       ENDIF

        write(lunmod,fmt='(1p,7e14.6)')
     >      F011,-f012, f033,-f034, phy,phz , Cstrn

      ENDIF

      RETURN

      END
