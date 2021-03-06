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
      SUBROUTINE FITSTA(IO,FITIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL FITIN, FITING
      SAVE FITING
      SAVE NUMKLE

      LOGICAL FITBYI, FITBYO, FITBYD
      SAVE FITBYD

      LOGICAL FITFNL, FITFNI, FITFNO
      SAVE FITFNL

      LOGICAL FITRBL, FITRBI, FITRBO
      SAVE FITRBL

      LOGICAL FITLST, FITLSI, FITLSO
      SAVE FITLST

      LOGICAL FITNHI, FITNHO, FITNHB

      DATA FITING / .FALSE. /
      DATA NUMKLE / 0 /
      DATA FITBYD / .FALSE. /
      DATA FITFNL / .FALSE. /
      DATA FITRBL / .FALSE. /
      DATA FITLST / .FALSE. /
      DATA FITNHB / .FALSE. /

C FITIN set to .true. (in zgoubi.f) when FIT procedure run
      IF(IO.EQ.5) THEN
C-------- read status
         FITIN = FITING
      ELSEIF(IO.EQ.6) THEN
C-------- save status
         FITING = FITIN
      ENDIF
      RETURN

      ENTRY FITST1(
     >             NUMKL)
         NUMKL = NUMKLE
      RETURN

      ENTRY FITST2(NUMKL)
         NUMKLE = NUMKL
      RETURN

! .T. if FIT has been completed, and pgm executing beyond keyword FIT[2}
      ENTRY FITST3(
     >             FITBYO)
         FITBYO = FITBYD
      RETURN
      ENTRY FITST4(FITBYI)
         FITBYD = FITBYI
      RETURN

! .T. if request for executing Zgoubi w new variables (following completion of FIT[2])
      ENTRY FITST5(
     >             FITFNO)
      FITFNO = FITFNL
      RETURN
      ENTRY FITST6(FITFNI)
      FITFNL = FITFNI
      RETURN

! .T. if FIT[2] embedded in REBELOTE
      ENTRY FITST7(
     >             FITRBO)
         FITRBO = FITRBL
      RETURN
      ENTRY FITST8(FITRBI)
         FITRBL = FITRBI
      RETURN

! .T. if FIT[2] is doing the last pass following completion (that can only happen if FITFNL=T)
      ENTRY FITST9(
     >             FITLSO)
         FITLSO = FITLST
      RETURN
      ENTRY FITSTX(FITLSI)
         FITLST = FITLSI
      RETURN

! .T. if FIT[2] is inhibited (by TWISS for instance, when FIT precedes FIT for orbit finding)
      ENTRY FITSTB(
     >             FITNHO)
         FITNHO = FITNHB
      RETURN
      ENTRY FITSTC(FITNHI)
         FITNHB = FITNHI
      RETURN

      END
