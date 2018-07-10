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
C  USA
C  -------
      SUBROUTINE RMCOBJ
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     *******************************************
C     READS DATA FOR OBJECT DEFINITION BY 'OBJET'
C     *******************************************
      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TXT

C----- BORO
      READ(NDAT,*) A(NOEL,1)
C----- KOBJ
      READ(NDAT,*) KOBJ
      A(NOEL,10) = KOBJ
C----- IMAX
      READ(NDAT,*) A(NOEL,20)
C----- Type of p.d.f, KY, KT, KZ, KP, KX, KD
      READ(NDAT,*) TXT
C----- Central values Yo, To, Zo, Po, Xo, Do
      READ(NDAT,*) TXT

      GOTO (1,2,3) KOBJ
      WRITE(NRES,FMT=
     >'(10X,'' ** SBR RMCOBJ:  no such object KOBJ='',I2)') KOBJ
      WRITE(NRES,FMT=
     >'(10X,''      KOBJ should be 1, 2 or 3'')')
      STOP
 
 1    CONTINUE
C----- dY, dT, dZ,...
      READ(NDAT,*) TXT
C----- Cut-offs Nsy, Nst, Nsz,...
      READ(NDAT,*) TXT
C----- No, Co, C1, C2 C3
      READ(NDAT,*) (A(NOEL,I),I=70,74)
C----- IR1, IR2, IR3
      READ(NDAT,*) (A(NOEL,I),I=80,82)
      RETURN
 
 2    CONTINUE
C----- # of bars IY, IT , IZ,...
      READ(NDAT,*) TXT
C----- Distances between bars PY, PT, PZ,...
      READ(NDAT,*) TXT
C----- Widths of the bars dY, dT, dZ,...
      READ(NDAT,*) TXT
C----- Cut-offs Nsy, Nst, Nsz,...
      READ(NDAT,*) TXT
C----- No, Co, C1, C2 C3
      READ(NDAT,*) (A(NOEL,I),I=90,94)
C----- IR1, IR2, IR3
      READ(NDAT,*) (A(NOEL,I),I=95,97)
      RETURN
 
 3    CONTINUE
      READ(NDAT,*) (A(NOEL,I),I=50,53)
      IF(A(NOEL,53) .LT. 0.D0) THEN
        BACKSPACE(NDAT,ERR=99)
        READ(NDAT,*) (A(NOEL,I),I=50,54)
        A(NOEL,53) = -A(NOEL,53) 
      ELSE
        A(NOEL,54) = 0.D0
      ENDIF 
      READ(NDAT,*) (A(NOEL,I),I=60,63)
      IF(A(NOEL,63) .LT. 0.D0) THEN
        BACKSPACE(NDAT,ERR=99)
        READ(NDAT,*) (A(NOEL,I),I=60,64)
        A(NOEL,63) = -A(NOEL,63) 
      ELSE
        A(NOEL,64) = 0.D0
      ENDIF 
      READ(NDAT,*) (A(NOEL,I),I=70,73)
      IF(A(NOEL,73) .LT. 0.D0) THEN
        BACKSPACE(NDAT,ERR=99)
        READ(NDAT,*) (A(NOEL,I),I=70,74)
        A(NOEL,73) = -A(NOEL,73) 
      ELSE
        A(NOEL,74) = 0.D0
      ENDIF 
      READ(NDAT,*) (A(NOEL,I),I=80,82)
      RETURN
 99   STOP ' *** Error upon BACKSPACE in SUBROUTINE RMCOBJ'
      END
