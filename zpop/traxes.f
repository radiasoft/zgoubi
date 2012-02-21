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
      SUBROUTINE TRAXES(XMI,XMA,YMI,YMA,MOD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE PROP
      DATA PROP /1.D0/

C----- CORR has been introduced so to make sure plot window does encompass xy-min/max
C         ... seems useless !
      CORR = 0.D0      !!!!!!(XMA-XMI)*1.D-1
      XM1 = XMI - CORR
      XM2 = XMA + CORR
      CORR = 0.D0     !!!!!!!!(YMA-YMI)*1.D-1
      YM1 = YMI - CORR
      YM2 = YMA + CORR
C        call fbgtxt
C        write(*,*) ' ************************* sbr traxes, PROP= ',prop
      IF(PROP.NE.1.D0) THEN
C------- Proportional axis (a circle should appear as a circle)
        DX = XM2-XM1
        DXP = DX*PROP
        XM1 = XM1-(DXP-DX)/2.D0
        XM2 = XM2+(DXP-DX)/2.D0
      ENDIF
C      XM1 = XMI
C      XM2 = XMA
C      IF(PROP.NE.1.D0) THEN
CC------- Proportional axis (a circle should appear as a circle)
C        DX = XM2-XM1
C        DXP = DX*PROP
C        XM1 = XM1-(DXP-DX)/2.D0
C        XM2 = XM2+(DXP-DX)/2.D0
C      ENDIF
      CALL TRAXE(SNGL(XM1),SNGL(XM2),SNGL(YM1),SNGL(YM2),MOD)
      CALL FBGTXT
      RETURN
      ENTRY TRAXPI
      PROP = (1024.D0/768.D0)
      RETURN
      ENTRY TRAXPJ
      PROP = 1.D0
      RETURN
      END
