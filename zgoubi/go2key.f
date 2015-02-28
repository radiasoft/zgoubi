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
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE GO2KEY(NUMKEY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************
C     Places pointer right before key # numkey in zgoubi.dat
C     ******************************************************
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      CHARACTER(40) TXT40
      INTEGER DEBSTR, FINSTR

      REWIND(NDAT)

      NK = 0
 1    CONTINUE
        READ(NDAT,fmt='(A)',END=98,ERR=99) TXT40
        IF(DEBSTR(TXT40) .GE. 1) THEN
          TXT40 = TXT40(DEBSTR(TXT40):FINSTR(TXT40))
        ELSE
          do i = 1, 40
            TXT40(i:i) = ' '
          enddo
        ENDIF
C        WRITE(*,fmt='(A,2i4,A40)') ' GO2KEY ::::::: ',NUMKEY, NK ,TXT40
        IF(TXT40(1:1).NE.'''') GOTO 1
        NK = NK+1
        IF(NK .LT. NUMKEY) GOTO 1

      BACKSPACE(NDAT)
C      WRITE(*,*) ' GO2KEY ::::::: ',NUMKEY, NK ,TXT40
 98   RETURN
 99   STOP ' *** ERROR IN SBR GO2KEY'
      END 
