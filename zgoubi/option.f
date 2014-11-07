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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory       
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE OPTION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER(80) TA
      PARAMETER (MXTA=45)
      COMMON/DONT/ TA(MXL,MXTA)
      CHARACTER(40) TXT1, TXT2
      INTEGER DEBSTR, FINSTR
      SAVE NRSAV
      SAVE KWROFF
      LOGICAL FITING
      
      DATA NRSAV / -11111 /
      data kwroff /  0 /

C Numb of options. NBOP lines should follow
      NY = NINT(A(NOEL,1))
      NBOP = NINT(A(NOEL,2))

      IF(NY*NBOP.EQ.0) THEN
        IF(NRSAV .EQ. -11111) THEN 
          IF(NRES.GT.0) WRITE(ABS(NRES),FMT='(/,T25,I2,A)') 
     >    ' ''OPTIONS''  is  inhibited,  thus  no  options  set.'
        ENDIF
        goto 99
      ENDIF

      IF(NBOP.GT.40)
     >CALL ENDJOB('SBR option : nmbr of options exceded ; max is ',40)
       
      IF(NRSAV .EQ. -11111) THEN 
          IF(NRES.GT.0) WRITE(ABS(NRES),FMT='(T25,I2,A)') 
     >    NBOP, ' option(s) expected.  Option(s) found :'
      ENDIF

      DO I = 1, NBOP
        READ(TA(NOEL,I),*) TXT1
        IF(NRES.GT.0) THEN
          IF(NRSAV .EQ. -11111) WRITE(ABS(NRES),FMT='(/,T5,A,I2,2A)') 
     >    'Option # ',I,' : ',
     >    TA(NOEL,I)(DEBSTR(TA(NOEL,I)):FINSTR(TA(NOEL,I)))
        ENDIF
        IF(TXT1(DEBSTR(TXT1):FINSTR(TXT1)) .EQ. 'WRITE') 
     >    READ(TA(NOEL,I),*) TXT1, TXT2
      ENDDO

      IF(NRSAV .EQ. -11111)  THEN 
         IF(NRES.GT.0) WRITE(ABS(NRES),FMT='(/,T25,A)') 
     >   ' Action(s) taken :' 
      ENDIF

      IF(TXT2(DEBSTR(TXT2):FINSTR(TXT2)) .EQ. 'OFF') THEN
 2      if(nrsav .eq. -11111) then
          IF(NRES.GT.0) THEN
            WRITE(ABS(NRES),FMT='(/,T5,A)') 
     >      'WRITE OFF -> '//
     >      'A lot of (almost all) WRITE statements will be inhibited !'
            WRITE(ABS(NRES),201)
 201        FORMAT(/,132('*'))
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
            WRITE(ABS(NRES),201)
            WRITE(ABS(NRES),*)  '  Keyword ''OPTIONS'' '
            WRITE(ABS(NRES),FMT='(T25,I2,A,/)') 
     >      NBOP, ' option(s) expected.  Option(s) found :'
            WRITE(ABS(NRES),FMT='(/,T5,A)') 'WRITE ON -> '//
     >      '''WRITE'' bit in ''OPTIONS'' set to 1.'
          ENDIF
          CALL REBEL1(
     >                 KWRT)
          IF(KWRT.NE.0)  KWROFF = 0
        ENDIF
      ENDIF

 99   RETURN

      ENTRY OPTIO1(
     >             KWROFO)
      KWROFO = KWROFF
      RETURN
      END
