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
      SUBROUTINE RSPN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ****************************
C     READS DATA FOR SPIN TRACKING
C     ****************************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
 
C     ... INITIAL SPIN DISTRIBUTION OPTION
      READ(NDAT,*) KSO
      A(NOEL,1) = KSO
 
      IF     (KSO .EQ. 4) THEN
        DO 1 I=1,IMAX
          READ(NDAT,*) (SI(J,I),J=1,3)
 1      CONTINUE
      ELSEIF(KSO .EQ. 5) THEN
C       ... TO, PO = MEAN INITIAL PRCESSION DIRECTION
        READ(NDAT,*) A(NOEL,10), A(NOEL,11)
C       ... AL, DA = CONE ANGLE AND D-ANGLE AROUND TO, PO
        READ(NDAT,*) A(NOEL,20), A(NOEL,21)
      ENDIF 
 
      RETURN
      END
