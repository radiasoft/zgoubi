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
      SUBROUTINE OPNMNL(LU2O,
     >                       NL,FN2O,OKOPN,CHANGE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKOPN, CHANGE
      CHARACTER FN2O*(*)

      COMMON/LUN/NDAT,NRES,NPLT,NFAI,NMAP,NSPN

      CHARACTER DEFN2O*20
      LOGICAL EMPTY, BINARY
      INTEGER DEBSTR, FINSTR
      CHARACTER(11) FRMT
      INCLUDE 'FILPLT.H'
      INCLUDE 'FILFAI.H'

      IF    (LU2O .EQ. NFAI) THEN
        DEFN2O = FILFAI
      ELSEIF(LU2O .EQ. NPLT) THEN
        DEFN2O = FILPLT
      ELSEIF(LU2O .EQ. NSPN) THEN
        DEFN2O = 'zgoubi.spn'
      ELSE
C        DEFN2O = 'none'
        IF( EMPTY(FN2O) ) THEN
          DEFN2O='none'
        ELSE
          DEFN2O = FN2O
        ENDIF
      ENDIF

      IF(OKOPN) THEN
        CLOSE(NL)
        OKOPN = .FALSE.        
      ENDIF
      CHANGE=.TRUE.
      NL = LU2O

      WRITE(6,FMT='(/,2A)')
     >' Give input data file name  - Default will be ',DEFN2O
      WRITE(6,FMT='(/,''  ("none" for exit)'')') 
      READ(5,FMT='(A)') FN2O
      IF( EMPTY(FN2O) ) THEN
        FN2O = DEFN2O
      ELSE
        FN2O=FN2O(DEBSTR(FN2O):FINSTR(FN2O))
      ENDIF
      IF(FN2O .NE. 'none') THEN
        WRITE(6,*)
        WRITE(6,*) ' Trying to open ',FN2O,' ...' 
        IDB = DEBSTR(FN2O)
        BINARY=FN2O(IDB:IDB+1).EQ.'B_' .OR. FN2O(IDB:IDB+1).EQ.'b_'
        IF(BINARY) THEN 
          FRMT='UNFORMATTED'
        ELSE
          FRMT='FORMATTED'
        ENDIF
        OPEN(UNIT=NL, FILE=FN2O,STATUS='OLD',IOSTAT=IOS,
     >            FORM=FRMT)
        IF(IOS .NE. 0) THEN
          WRITE(6,*) '  Cannot open ',FN2O
          FN2O = 'none'
        ELSE
          OKOPN = .TRUE.
          WRITE(6,FMT='(2A)') ' Opened file ',FN2O
          WRITE(6,*) ' Particle coordinates will be read from ',FN2O
          IDB=DEBSTR(FN2O)
          CALL HEADER(NL,4,BINARY,*49)
        ENDIF
        IF(LU2O.EQ.NFAI .OR. LU2O.EQ.NPLT .OR. LU2O.EQ.NSPN) 
     >                                          CALL READC8(BINARY)
      ELSE
        WRITE(6,*)
        WRITE(6,*) ' No  file  opened...'
      ENDIF
 49   RETURN
      END
