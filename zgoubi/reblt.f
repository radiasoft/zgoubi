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
      SUBROUTINE REBLT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ***************************************
C     READS DATA FOR FIT PROCEDURE WITH 'FIT'
C     ***************************************
      INCLUDE 'MXLD.H'
      SAVE NLA, NLB
      PARAMETER (KSIZ=10)
      SAVE KLERB

      DATA NLA, NLB / 1, MXL /
      DATA KLERB/ 11 /

      ENTRY REBLT4(KLERBI)
      KLERB = KLERBI
      RETURN

      ENTRY REBLT5(
     >             NLAO,NLBO)
      NLAO = NLA
      NLBO = NLB
      RETURN

      ENTRY REBLT6(NLAI,NLBI)
      NLA = NLAI
      NLB = NLBI
      RETURN

      ENTRY REBLT7(
     >             NLBO,KLERBO)
      NLBO = NLB
      KLERBO = KLERB
      RETURN

      END
