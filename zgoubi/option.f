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
C  -------
      SUBROUTINE OPTION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      PARAMETER (LNTA=132) ; CHARACTER(LNTA) TA
C      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      CHARACTER(40) TXT1, TXT2
      INTEGER DEBSTR, FINSTR
      SAVE NRSAV
      SAVE KWROFF
      LOGICAL FITING
      LOGICAL CONSTY
      
      DATA NRSAV / -11111 /
      DATA KWROFF /  0 /
      DATA CONSTY / .FALSE. /

C Numb of options. NBOP lines should follow
      NY = NINT(A(NOEL,1))
      NBOP = NINT(A(NOEL,2))

      IF(NY*NBOP.EQ.0) THEN
        IF(NRSAV .EQ. -11111) THEN 
          IF(NRES.GT.0) WRITE(ABS(NRES),FMT='(/,T25,A)') 
     >    ' ''OPTIONS''  is  inhibited,  no  option  will  be  set.'
        ENDIF
        GOTO 99
      ENDIF

      IF(NBOP.GT.40)
     >CALL ENDJOB('SBR option : nmbr of options exceded ; max is ',40)
       
      IF(NRSAV .EQ. -11111) THEN 
               IF(NRES.GT.0) WRITE(ABS(NRES),FMT='(T10,A,I0,A,/)') 
     >         'A list of ',NBOP,' option(s) is expected.  '//
     >         'List and actions taken are as follows :'
      ENDIF

C      DO I = 1, NBOP
C        READ(TA(NOEL,I),*,ERR=88,END=88) TXT1
C        IF(NRES.GT.0) THEN
C          IF(NRSAV .EQ. -11111) WRITE(ABS(NRES),FMT='(/,T5,A,I2,2A)') 
C     >    'Option # ',I,' : ',
C     >    TA(NOEL,I)(DEBSTR(TA(NOEL,I)):FINSTR(TA(NOEL,I)))
C        ENDIF
C        IF(TXT1(DEBSTR(TXT1):FINSTR(TXT1)) .EQ. 'WRITE') 
C     >    READ(TA(NOEL,I),*) TXT1, TXT2
C      ENDDO

      DO I = 1, NBOP

       READ(TA(NOEL,I),*,ERR=88,END=88) TXT1

       IF(NRES.GT.0) THEN
          IF(NRSAV .EQ. -11111) WRITE(ABS(NRES),FMT='(/,T5,A,I2,2A)') 
     >    'Option # ',I,' : ',
     >    TA(NOEL,I)(DEBSTR(TA(NOEL,I)):FINSTR(TA(NOEL,I)))
       ENDIF

       IF(TXT1(DEBSTR(TXT1):FINSTR(TXT1)) .EQ. 'WRITE') THEN

         READ(TA(NOEL,I),*) TXT1, TXT2

         IF(TXT2(DEBSTR(TXT2):FINSTR(TXT2)) .EQ. 'OFF') THEN
           IF(NRSAV .EQ. -11111) THEN
             IF(NRES.GT.0) THEN
               WRITE(ABS(NRES),FMT='(/,T5,A)') 
     >         'WRITE OFF -> A lot of (almost all) '//
     >         'WRITE statements will be inhibited !'
               WRITE(ABS(NRES),FMT='(/,132(''*''))')
             ENDIF
           ENDIF
           KWROFF = 1
           NRSAV = NRES
           NRES = -ABS(NRES)
         ELSE
           CALL FITSTA(5,
     >                FITING)
           IF(NRSAV .NE. -11111 .AND. .NOT. FITING) THEN 
C          IF(NRES.GT.0) THEN
C Yann, 14-03-07. Necessary for the online model to work
             IF(ABS(NRES).GT.0) THEN
               NRES = ABS(NRES)
               WRITE(ABS(NRES),FMT='(/,T5,A)') 'WRITE ON -> '//
     >         '''WRITE'' bit in ''OPTIONS'' set to 1.'
             ENDIF
             CALL REBEL1(
     >                 KWRT)
             IF(KWRT.NE.0)  KWROFF = 0
           ENDIF
         ENDIF

       ELSEIF(TXT1(DEBSTR(TXT1):FINSTR(TXT1)) .EQ. 'CONSTY') THEN

         READ(TA(NOEL,I),*) TXT1, TXT2
         
         CONSTY = TXT2(DEBSTR(TXT2):FINSTR(TXT2)) .EQ. 'ON'
         CALL INTEGA(CONSTY)
         CALL TRANS2(CONSTY)

         IF(ABS(NRES).GT.0) THEN
           NRES = ABS(NRES)
           WRITE(ABS(NRES),FMT='(5X,T5,A)') '- rays will be forced to'
     >     //' constant Y and Z in pgm integr.'
           WRITE(ABS(NRES),FMT='(5X,T5,A)') '- mid-plane symmetry'
     >     //' test in pgm transf (i.e., ''dejaca'' procedure) '
     >     //'is inhibited.'
           WRITE(ABS(NRES),FMT='(5X,T5,A)') '- coordinates and field '
     >     //'may be checked using IL = 1, 2 OR 7 '
c     >     //'will stored in file zgoubi.consty.out '
         ENDIF

       ENDIF
      ENDDO

 99   RETURN

      ENTRY OPTIO1(
     >             KWROFO)
      KWROFO = KWROFF
 88   RETURN
      END
