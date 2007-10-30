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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE SAVPLT(*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 FNAM
      INTEGER DEBSTR, FINSTR
      LOGICAL EMPTY, IDLUNI
 
 1    CONTINUE
      WRITE(6,*) '  GIVE FILE NAME  (DEFAULT WILL BE screen-save.ps) :'
      READ(5,100,ERR=1) FNAM
 100  FORMAT(A80)

      IF( EMPTY(FNAM) ) THEN
        FNAM = 'screen-save.ps'
      ELSE
        FNAM = FNAM( DEBSTR(FNAM) : FINSTR(FNAM) )
      ENDIF

      IF (IDLUNI(IUN)) THEN
C        WRITE(6,*) '  Logical unit ',IUN,' found free...',
C     >     ' Trying to connect to ',FNAM
        OPEN(UNIT=IUN,FILE=FNAM,STATUS='NEW',ERR=997)
      ELSE
        WRITE(6,*) ' *** Problem : No idle unit number ! '
        GOTO 998
      ENDIF

C      CALL HPLCAP(21)      ! Cern lib
      CALL SAVECR(FNAM)    ! lnslib

      CLOSE(IUN)
      WRITE(6,*) '  Done. Screen saved in ', FNAM
      RETURN

 997  WRITE(6,*) ' Error upon OPEN statement '
 998  WRITE(6,*) ' ***  NO  SCREEN SAVE...'

      RETURN 1
      END
