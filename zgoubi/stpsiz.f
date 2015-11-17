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
      SUBROUTINE STPSIZ(NDAT,NOEL,ND,
     >                               A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'MXLD.H'
      DIMENSION A(MXL,*)
      CHARACTER TXT*80, TXTEMP*80
      CHARACTER STRA(1)*80
      LOGICAL STRCON
      INTEGER DEBSTR, FINSTR

      DIMENSION KPSTO(MXL), IPRSTO(MXL)
      SAVE KPSTO, IPRSTO

      READ(NDAT,FMT='(A)') TXT
      TXT = TXT(DEBSTR(TXT):FINSTR(TXT))
      CALL STRGET(TXT,1,
     >                     IDUM,STRA)
      TXT = STRA(1)
      IFIN = FINSTR(TXT)
      KPAS = 0
      IPREC = 0

      IF(STRCON(TXT,'E',IS) .OR. STRCON(TXT,'e',IS)) THEN
        READ(TXT,*) TEMP
        IF(TEMP.GT.1.D10) THEN
C--------- Old method for entrance|body|exit number of steps -> convert to new method
          TXT = TXT(1:IS-1)
          IF(STRCON(TXT,'.',IS)) WRITE(TXTEMP,FMT='(A)') 
     >       '#'//TXT(IS+1:IS+3)//'|'//TXT(1:IS-1)//'|'//TXT(IS+1:IS+3)
          TXT = TXTEMP(DEBSTR(TXTEMP):FINSTR(TXTEMP))
          IFIN = FINSTR(TXT)
        ENDIF
      ELSEIF(STRCON(TXT,'V',IS) .OR. STRCON(TXT,'v',IS)) THEN
C------- Variable number of steps in entrance|body|exit regions
        KPAS = 2
        READ(TXT(IS+2:IS+3),*) IPREC
        TXT = TXT(1:IS-1)
        IF(STRCON(TXT,'.',IS)) WRITE(TXTEMP,FMT='(A)') 
     >       '#'//TXT(IS+1:IS+3)//'|'//TXT(1:IS-1)//'|'//TXT(IS+1:IS+3)
        TXT = TXTEMP(DEBSTR(TXTEMP):FINSTR(TXTEMP))
        IFIN = FINSTR(TXT)
      ENDIF

      IF(TXT(1:1) .EQ. '#') THEN
C------- Coded step size with form  #stp_1|stp_2|stp_3, stp_i= arbitrary integer
C                                    entr body exit
        IF(IFIN.LT.6) GOTO 99
        IF(STRCON(TXT,'|',IS)) THEN
          IF(IS.GT.IFIN-3) GOTO 99
          READ(TXT(2:IS-1),*) A(NOEL,ND)
          IS1 = IS
          IF(STRCON(TXT(IS1+1:IFIN),'|',IS)) THEN
            IF(IS1+IS.GT.IFIN-1) GOTO 99
            READ(TXT(IS1+1:IS1+IS-1),*) A(NOEL,ND+1)
            READ(TXT(IS1+IS+1:IFIN),*) A(NOEL,ND+2)
          ENDIF            
        ELSE
          GOTO 99
        ENDIF
C Modif, FM, Dec. 05
C        ND = ND+2
        KPAS = 1
        IPREC = 0
      ELSE
        READ(TXT,*) A(NOEL,ND)
      ENDIF
      
      KPSTO(NOEL) = KPAS
      IPRSTO(NOEL) = IPREC
      CALL CHXC1W(KPAS,IPREC)
      RETURN


      ENTRY STPSI1(NUML)
      CALL CHXC1W(KPSTO(NUML),IPRSTO(NUML))
      RETURN

 99   CONTINUE
      STOP ' *** DATA ERROR : Upon reading coded step size ***  '

      END
