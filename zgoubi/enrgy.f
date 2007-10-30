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
      SUBROUTINE ENRGY(ERROR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      COMMON/PTICUL/ AM,Q,G,TO
      COMMON/RIGID/ BORO,DPREF,DP,BR
      SAVE R, V
      SAVE UX
C----- Momentum is in MeV/c (BR=BORO*DP in kG.cm, see SBR MAJTRA)
C      P = BR *CL*1.D-9*Q/QE
C      P = BR *CL*1.D-9*Q
      P = BR *CL9*Q
      E = P*P/(2.D0*AM) + V
C      P0 = BORO*CL*1.D-9*Q/QE
C      P0 = BORO*CL*1.D-9*Q
      P0 = BORO*CL9*Q
      E0 = P0*P0/(2.D0*AM) 
      ERROR = (E - E0)/E0
      RVX = R * UX* P
      RETURN
      ENTRY ENRGW(RIN,VIN)
      R = RIN
      V = VIN
      RETURN
      ENTRY VTHET(UXIN)
      UX = UXIN
      RETURN
      END
