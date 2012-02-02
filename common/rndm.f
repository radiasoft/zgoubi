C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory   
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      FUNCTION RNDM(IR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      CHARACTER TXT*16
c      INTEGER DEBSTR
c      REAL R
      logical first 
      save first
C year month day time-utc hour min secon ms
      dimension ival(8)
      data first / .true. /

      if(first) then
        first = .false.
        call date_and_time(VALUES=ival)
        ir = (1+ival(5))*100000 + ival(6)*1000 + ival(7)*10 +ival(8)
        rndm = rand(ir) 
c        call srand(ir)
c        write(*,*) ' rndm ',ival(5),ival(6),ival(7),ival(8),ir,rndm
      endif
      rndm = rand(0) 

c      CALL RANDOM_NUMBER(R)
c      RNDM=R
c      WRITE(TXT,FMT='(E14.6)') RNDM
c      READ(TXT(DEBSTR(TXT)+2:DEBSTR(TXT)+8),FMT='(I6)') IR      
      
      RETURN
      END
