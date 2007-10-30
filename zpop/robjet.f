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
      SUBROUTINE ROBJET
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     *******************************************
C     READS DATA FOR OBJECT DEFINITION BY 'OBJET'
C     *******************************************
      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,20)
      PARAMETER (MXJ=7)
      DIMENSION IA(5)

      PARAMETER(MXJ1=MXJ-1)

C----- BORO
      READ(NDAT,*,ERR=99) A(NOEL,1)
C----- KOBJ
      READ(NDAT,*,ERR=99) A(NOEL,10)
 
      K=A(NOEL,10)
      IF(K .LT. 0) K=-K

      GOTO (1,2,3,1,5,5,7,8) K

      WRITE(NRES,FMT=
     >'(10X,''SBR ROBJET:    NO  SUCH  OBJECT  KOBJ='',I2)') K
      STOP
 
 1    CONTINUE
      READ(NDAT,*,ERR=99) (IA(I),I=1,MXJ1)
      DO 11 I=1,MXJ1
 11     A(NOEL,19+I) = IA(I)
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=30,35)
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=40,45)
      RETURN
 
 2    CONTINUE
      READ(NDAT,*,ERR=99) (IA(I),I=1,2)
      A(NOEL,20) = IA(1)
      A(NOEL,21) = IA(2)
      DO 21 I=1,IA(1)
        READ(NDAT,*,ERR=99)  DUM
 21   CONTINUE
      READ(NDAT,*,ERR=99) (IEX,I=1,IA(1))
      RETURN
 
 3    CONTINUE
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=20,21)
      READ(NDAT,100,ERR=99) TA(NOEL,1)
 100  FORMAT(A80)
      RETURN
 
 5    CONTINUE
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=20,25)
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=30,35)
      RETURN
 
 7    CONTINUE
      RETURN
 
 8    CONTINUE
C----- IY, IZ, IX
      READ(NDAT,*,ERR=99) (A(NOEL,19+I),I=1,3)
C----- Center of ellipsoid (Y, T, Z, P, X, D
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=30,35)
C----- alpha, beta, epsilon/pi for Y, Z, X phase-spaces
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=40,42)
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=50,52)
      READ(NDAT,*,ERR=99) (A(NOEL,I),I=60,62)
      RETURN
 
 99   WRITE(6,*) 
     >  ' *** Execution stopped upon READ : invalid input in OBJET'
      WRITE(NRES ,*) 
     >  ' *** Execution stopped upon READ : invalid input in OBJET'
      STOP
      END
