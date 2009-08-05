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
      SUBROUTINE SCUMUL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE "MAXCOO.H"
      COMMON/UNITS/ UNIT(MXJ)

      DATA XL, SCUM, TCUM / 0.D0, 0.D0, 0.D0 / 
      SAVE XL, SCUM, TCUM

      ENTRY SCUMW(XLI)
      XL=XLI
      SCUM = SCUM + XL
C----- Compute cumulative time. Default is for proton. 
      AAM = AM
      QQ = Q
      IF(AAM .EQ. 0.D0) AAM = AMPROT
      IF(QQ .EQ. 0.D0) QQ = QE
C      PREF = (BORO*DPREF) *CL*1.D-9*QQ/QE
C      PREF = (BORO*DPREF) *CL*1.D-9*QQ
      PREF = (BORO*DPREF) *CL9*QQ
      BTA = PREF / SQRT(PREF*PREF+AAM*AAM)
C----- XL is in centimeters, TCUM is seconds
      DT = XL / (CL*BTA) *UNIT(5)
      TCUM = TCUM + DT
      RETURN

      ENTRY SCUMR(
     >            XLO,SCUMO,TCUMO)
      XLO=XL
      SCUMO = SCUM
      TCUMO = TCUM
      RETURN

      ENTRY SCUMS(SETS)
      SCUM = SETS
C      TCUM = SETT
      RETURN
      END
