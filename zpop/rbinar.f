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
      SUBROUTINE RBINAR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     *********************
C     READS DATA FOR BINARY
C     *********************
      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      PARAMETER (MXTA=20)
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,MXTA)
 
C     ... # OF FILES TO TRANSLATE
      READ(NDAT,*) NF
      IF(NF.GT.MXTA) STOP ' IN SBR RBINAR: TO MANY FILES'
      A(NOEL,1)=NF
C     ... FILE NAMES
      DO 1 I=1,NF
        READ(NDAT,200) TA(NOEL,I)
 200    FORMAT(A)
 1    CONTINUE
 
      RETURN
      END
