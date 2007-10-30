C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <meot@lpsc.in2p3.fr>
C  Service Acc�lerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE SPNROT(I,AX,AY,AZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------------------------------------------------------
C     TOURNE LE REFERENTIEL DU SPIN PAS APRES
C     PAS, DANS CHANGREF ET DANS LES AIMANTS EN COORD. POLAIRES.
C     ---------------------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
 
      IF(AX .NE. 0.D0) THEN
        CA = COS(AX)
        SA = SIN(AX)
        SF2 = SF(2,I)
        SF(2,I) =  CA*SF2 + SA*SF(3,I)
        SF(3,I) = -SA*SF2 + CA*SF(3,I)
      ENDIF
      IF(AY .NE. 0.D0) THEN
        CA = COS(AY)
        SA = SIN(AY)
        SF3 = SF(3,I)
        SF(3,I) =  CA*SF3 + SA*SF(1,I)
        SF(1,I) = -SA*SF3 + CA*SF(1,I)
      ENDIF
      IF(AZ .NE. 0.D0) THEN
        CA = COS(AZ)
        SA = SIN(AZ)
        SF1 = SF(1,I)
        SF(1,I) =  CA*SF1 + SA*SF(2,I)
        SF(2,I) = -SA*SF1 + CA*SF(2,I)
      ENDIF

      RETURN
      END
