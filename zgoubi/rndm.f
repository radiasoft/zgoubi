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
      FUNCTION RNDM()
CDescription:
C    RAND(FLAG) returns a pseudo-random number from a uniform distribution between 0 and 1. 
C  If FLAG is 0, the next number in the current sequence is returned; if FLAG is 1, the generator 
C  is restarted by CALL SRAND(0); if FLAG has any other value, it is used as a new seed with SRAND.
C    This intrinsic routine is provided for backwards compatibility with GNU Fortran 77. It implements 
C  a simple modulo generator as provided by g77. For new code, one should consider the use of 
C    RANDOM_NUMBER as it implements a superior algorithm.
CStandard:
C    GNU extension
Class:
C    Function
CSyntax:
C    RESULT = RAND(I)
CArguments:
C    I 	Shall be a scalar INTEGER of kind 4.
CReturn value:
C    The return value is of REAL type and the default kind.
CExample:
C          program test_rand
C            integer,parameter :: seed = 86456         
C            call srand(seed)
C            print *, rand(), rand(), rand(), rand()
C            print *, rand(seed), rand(), rand(), rand()
C          end program test_rand

C For ifort compiler. Comment otherwise
c        use ifport

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      CHARACTER TXT*16
c      INTEGER DEBSTR
      REAL R
      logical first 
      save first
C year month day time-utc hour min secon ms
c      dimension ival(8)
c      save ir
      dimension i1(12)

      data first / .true. /
      data ir / 123456 /
      data i1 / 12 * 123456 / 

C For ifort compiler. Can be commented otherwise. 

      if(first) then
        first = .false.
C--------------------------
C Chose here between random seed or not
c        call date_and_time(VALUES=ival)
c        ir = (1+ival(5))*100000 + ival(6)*1000 + ival(7)*10 +ival(8)

c        call srand(ir)
c        CALL init_random_seed()

C--------------------------
c        write(*,*) ' rndm ',ival(5),ival(6),ival(7),ival(8),ir,rndm
      endif

c      R = rand()
      CALL RANDOM_NUMBER(R)

      RNDM=R
 
c          write(*,*) ' rndm ',r,ir,rndm
c             read(*,*)

c      WRITE(TXT,FMT='(E14.6)') RNDM
c      READ(TXT(DEBSTR(TXT)+2:DEBSTR(TXT)+8),FMT='(I6)') IR      
      
      RETURN

      entry rndm2(iri)
C      call srand(iri)
       i1(1) = iri
        call RANDOM_seed(put=i1)
      return

      END
