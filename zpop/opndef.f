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
      SUBROUTINE OPNDEF(LU2O,DEFN2O,LUO,
     >                                  FNO,OKOPN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) DEFN2O, FNO
      LOGICAL OKOPN

      LOGICAL EXS, OPN, BINARY 
      CHARACTER(11) FRMT
      INCLUDE "FILHDF.H"
      
      IF( LU2O .EQ. -1) THEN
C------- Looks for a free LUO starting from #10
        OKOPN = .FALSE.
        LU2O = 9
 1      CONTINUE  
        LU2O = LU2O+1
        IF(LU2O.EQ.100) GOTO 96
        INQUIRE(UNIT=LU2O,EXIST=EXS,OPENED=OPN,IOSTAT=IOS)
        IF( IOS .EQ. 0 .AND. OPN ) GOTO 1
      ELSEIF(LU2O .EQ. LUO) THEN
C------- Check whether LUO is already open
        INQUIRE(UNIT=LUO,EXIST=EXS,OPENED=OPN,NAME=FNO,IOSTAT=IOS)
        IF(IOS .EQ. 0) THEN
          OKOPN = OPN
        ELSE
          OKOPN = .FALSE.
        ENDIF
      ELSE
        IF(LUO.GT.0) CLOSE(LUO)
        OKOPN = .FALSE.
      ENDIF

      IF(.NOT. OKOPN) THEN
C--------- Check existence of DEFN2O
        INQUIRE(FILE=DEFN2O,EXIST=EXS,IOSTAT=IOS)

        IF(IOS .EQ. 0) THEN

          BINARY=DEFN2O(1:2).EQ.'B_' .OR. DEFN2O(1:2).EQ.'b_'
          IF(BINARY) THEN 
            FRMT='UNFORMATTED'
          ELSE
            FRMT='FORMATTED'
          ENDIF
          IF(EXS) THEN
            OPEN(UNIT=LU2O,FILE=DEFN2O,STATUS='OLD',ERR=99,IOSTAT=IOS,
     >           FORM=FRMT)
            IF(IOS.NE.0) GOTO 97
            CALL HEADER(LU2O,4,BINARY,*99)
          ELSE
            OPEN(UNIT=LU2O,FILE=DEFN2O,STATUS='NEW',ERR=99,IOSTAT=IOS,
     >           FORM=FRMT)
            IF(IOS.NE.0) GOTO 97
          ENDIF

          FNO = DEFN2O
          LUO = LU2O
          OKOPN = .TRUE.
        ELSE

          WRITE(6,*) ' Exec error occured in Subroutine OPNDEF : '
          WRITE(6,*) ' IOS NON-ZERO '

        ENDIF
C        IF(OKOPN) WRITE(6,FMT='(2A)') 'Opened file is ',DEFN2O
C      ELSE
C         WRITE(6,*) ' Already opened file ',DEFN2O
        CALL READC8(BINARY)
      ENDIF

      IF(OKOPN) WRITE(6,FMT='(2A)') 'Opened file is ',DEFN2O

      RETURN

 96   WRITE(6,*) ' Exec error occured in Subroutine OPNDEF : ' 
      WRITE(6,*) ' Logical unit # exceeds 100 ; CANNOT OPEN ',DEFN2O
      RETURN

 97   WRITE(6,*) ' Error occured while in Subroutine OPNDEF : '
      WRITE(6,*) '      IOS NON-ZERO ; CANNOT OPEN ',DEFN2O
 99   CONTINUE
      RETURN

C 98   CONTINUE
C      CLOSE(LU2O)
C      RETURN

      END
