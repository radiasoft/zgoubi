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
      SUBROUTINE ET2RES(LUNW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IDLUNI, OK, GTTEXT, EMPTY
      INTEGER DEBSTR, FINSTR
      CHARACTER(200) TXT200

      SAVE ALF1, ALF2, BET1, BET2, PHY, PHZ, CSTRN

      OK = IDLUNI(
     >            NL)
      IF(OK) THEN
        OPEN(UNIT=NL,FILE='ETparam.res',ERR=91)
      ELSE
        GOTO 91
      ENDIF

C      WRITE(NRES,FMT='(A)') 

      M1 = -1
      OK = GTTEXT(M1,NL,'###########',
     >                              TXT200)

      OK = GTTEXT(M1,NL,'FRACTIONAL PART',
     >                              TXT200)
      READ(TXT200(70: 89),*,ERR=97,END=97) PHY
      READ(TXT200(90:106),*,ERR=97,END=97) PHZ

      OK = GTTEXT(M1,NL,'EDWARDS-TENG',
     >                              TXT200)      
      READ(NL,FMT='(A)',ERR=98,END=98) TXT200
      READ(TXT200(70: 89),*,ERR=97,END=97) ALF1
      READ(TXT200(90:106),*,ERR=97,END=97) ALF2
      READ(NL,FMT='(A)',ERR=98,END=98) TXT200
      READ(TXT200(70: 89),*,ERR=97,END=97) BET1
      READ(TXT200(90:106),*,ERR=97,END=97) BET2

      OK = GTTEXT(M1,NL,'- COUPLING STRENG',
     >                              TXT200)      
      READ(TXT200(85:99),*,ERR=97,END=97) CSTRN

      GOTO 99

 91   CONTINUE
      WRITE(NRES,*)      
     >'      ' //
     >' SBR ET2RES. ' // 
     >'COULD NOT OPEN FILE "ETparam.res", PROCEEDING WITHOUT...'
      PHY = 0.D0
      PHZ = 0.D0
      ALF1 = 0.D0
      ALF2 = 0.D0
      BET1 = 0.D0
      BET2 = 0.D0
      GOTO 99

 98   CONTINUE
      WRITE(NRES,*)      
      WRITE(NRES,*)      
     >'      ' //
     >' SBR ET2RES. FINISHED COPYING FROM ETparam.res TO ZGOUBI.RES ' //
     >'UPON EOR OR EOF. '
      PHY = 0.D0
      PHZ = 0.D0
      ALF1 = 0.D0
      ALF2 = 0.D0
      BET1 = 0.D0
      BET2 = 0.D0
      GOTO 99

 97   CONTINUE
      WRITE(NRES,*)      
      WRITE(NRES,*)      
     >'      ' //
     >' SBR ET2RES. EOR UPON BAD DATA IN TXT200. ' 
      PHY = 0.D0
      PHZ = 0.D0
      ALF1 = 0.D0
      ALF2 = 0.D0
      BET1 = 0.D0
      BET2 = 0.D0
      GOTO 99

 99   CONTINUE
      IF(OK) CLOSE(NL)
      RETURN

      ENTRY ET2RE1(
     >             F011,F012,F033,F034,PHYO,PHZO,CSTRNO)
      F011 = BET1
      F012 = -ALF1
      F033 = BET2
      F034 = -ALF2
      PHYO = PHY
      PHZO = PHZ
      CSTRNO = CSTRN
      RETURN

      END
