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
      SUBROUTINE RSOLEN(NDAT,NOEL,MXL,A,ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MXL,*)
C     -----------------------
C     READS DATA FOR SOLENOID
C     -----------------------
      INCLUDE "C.MARK.H"     ! COMMON/MARK/ KART,KALC,KERK,KUASEX

      CHARACTER(132) TXT132
      CHARACTER(30) XRBM(4)
      INTEGER DEBSTR
      LOGICAL STRCON
       
C----- IL
      READ(NDAT,*) A(NOEL,1)

C----- XL, RO, BO [, MODL]. Default is MODL=1
C      READ(NDAT,*) A(NOEL,10),A(NOEL,11),A(NOEL,12)
      READ(NDAT,FMT='(A)') TXT132
      IF(STRCON(TXT132,'!',
     >                     II))
     >TXT132 = TXT132(DEBSTR(TXT132):II-1)
      I4 = 4
      CALL STRGET(TXT132,I4,
     >                      I34,XRBM)
      READ(XRBM(1),*,ERR=98,END=98) A(NOEL,10)
      READ(XRBM(2),*,ERR=98,END=98) A(NOEL,11)
      READ(XRBM(3),*,ERR=98,END=98) A(NOEL,12)
      IF    (I34.EQ.3) THEN 
         A(NOEL,13) = 1
      ELSEIF(I34.EQ.4) THEN 
         READ(XRBM(4),*) A(NOEL,13)
      ELSE
         GOTO 98
      ENDIF
C----- XE, XS
      READ(NDAT,*) A(NOEL,20),A(NOEL,21)

C----- Integr. step
      ND = 30
      CALL STPSIZ(NDAT,NOEL,ND,
     >                         A)

      READ(NDAT,*) IA,(A(NOEL,I),I=ND+10+1,ND+10+3)
      A(NOEL,ND+10) = IA
 
      RETURN

 98   CONTINUE
      CALL ENDJOB('SBR RSOLEN. Input data error at line #',2)
      RETURN
      END
