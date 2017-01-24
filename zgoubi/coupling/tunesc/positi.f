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
C     written by Frédéric Desforges, frederic.desforges@grenoble-inp.org
C
      SUBROUTINE POSITI(NRES,ARCLEN,LABEL,KEYWOR)

      IMPLICIT DOUBLE PRECISION (A-H,M-Z)
      INTEGER NRES,IS  !,I,DEBSTR,FINSTR
      CHARACTER(300) BUFFER,LABEL,KEYWOR
      LOGICAL STRCON
 
 1    CONTINUE
      
      READ(NRES,FMT='(a)',ERR=96,END=97) BUFFER
      
      IF(STRCON(BUFFER,'Reference particle (#',IS) .
     >EQV. .TRUE.) THEN
         
c         WRITE(*,*) ' positi.  ', BUFFER
cc             read(*,*)

         READ(BUFFER(44:54),*) ARCLEN
         
 2       CONTINUE
         
         READ(NRES,FMT='(a)',ERR=96,END=97) BUFFER

         IF(STRCON(BUFFER,'*********************************************
     >******************************************************************
     >*****************',IS) .EQV. .TRUE.) THEN
            READ(NRES,FMT='(a)',ERR=96,END=97) BUFFER
            READ(BUFFER(30:39),*,ERR=96,END=97) KEYWOR
            READ(BUFFER(40:60),*,ERR=96,END=97) LABEL
c            WRITE(*,*) ' positi ',arclen,KEYWOR(1:9),' ',LABEL(1:20)
cc                 read(*,*)
            GOTO 97

         ENDIF

         GOTO 2
     
      ENDIF

      GOTO 1

 96   CONTINUE

      WRITE(*,FMT='(/,/,''ERROR IN READING OF ZGOUBI.RES'',/,/)')
      WRITE(*,*) 'PBM WITH ',KEYWOR(1:9),' ',LABEL(1:20)

 97   CONTINUE
      
      RETURN

       
      END
