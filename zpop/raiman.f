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
      SUBROUTINE RAIMAN(NDAT,NOEL,MXL,A,ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     *****************************
C     READS DATA FOR AIMANT, DIPOLE
C     *****************************
      DIMENSION A(MXL,*)
 
      READ(NDAT,*) (A(NOEL,I),I=1,3)
      READ(NDAT,*) (A(NOEL,I),I=4,5)
      READ(NDAT,*) (A(NOEL,I),I=6,9)
      READ(NDAT,*) (A(NOEL,I),I=10,14)
C       ... FACE ENTREE
      READ(NDAT,*) (A(NOEL,I),I=15,16)
      READ(NDAT,*) IIA,(A(NOEL,I),I=18,24)
      A(NOEL,17) = IIA
      READ(NDAT,*) (A(NOEL,I),I=25,30)
C       ... FACE SORTIE
      READ(NDAT,*) (A(NOEL,I),I=31,32)
      READ(NDAT,*) IIA,(A(NOEL,I),I=34,40)
      A(NOEL,33) = IIA
      READ(NDAT,*) (A(NOEL,I),I=41,46)
C       ... FACE LATERALE
      IF(A(NOEL,1) .EQ. 3) THEN
        READ(NDAT,*) (A(NOEL,I),I=47,48)
        READ(NDAT,*) IIA,(A(NOEL,I),I=50,56)
        A(NOEL,49) = IIA
        READ(NDAT,*) (A(NOEL,I),I=57,63)
      ENDIF
      READ(NDAT,*) IIA
      A(NOEL,64) = IIA
C      IF(IIA .NE. 0) STOP ' NO SHIMS AVAILABLE: TO BE IMPLEMENTED'
 
C     ... IRD
      READ(NDAT,*) A(NOEL,93)
C     ... XPAS
      READ(NDAT,*) A(NOEL,94)
      ND=94
C     ... KP
      READ(NDAT,*) IIA
      A(NOEL,95) = IIA
      IF( IIA .EQ. 2 ) THEN
        READ(NDAT,*) (A(NOEL,I),I=96,99)
      ELSE
        READ(NDAT,*) A(NOEL,96)
      ENDIF
 
      RETURN
      END
