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
      SUBROUTINE RREBEL(LABEL,KLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ***************************************
C     READS DATA FOR FIT PROCEDURE WITH 'FIT'
C     ***************************************
      CHARACTER(*) KLE(*)
      INCLUDE 'MXLD.H'
      CHARACTER(*) LABEL(MXL,2)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL

      CHARACTER(8) TXTA, TXTB
      PARAMETER(I300=300)
      CHARACTER(I300) TXT300, TXTMP
      INTEGER DEBSTR, FINSTR
      LOGICAL STRCON
      PARAMETER (I2=2, I4=4)
      CHARACTER(I300) STRA(I4)
      CHARACTER(30) STRING
      LOGICAL OKKLE

      PARAMETER (MXPRM=10, MXLST=4000)
      DIMENSION PARAM(MXPRM,MXLST)
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) TPRM(MXPRM,3)
      LOGICAL ISNUM
      PARAMETER (I3=3)

      DATA IA4 / 0 /

      LINE = 1
      READ(NDAT,FMT='(A)',ERR=98,END=98) TXT300
      IF(STRCON(TXT300,'!',
     >                     II))
     >  TXT300 = TXT300(DEBSTR(TXT300):II-1)
      CALL STRGET(TXT300,I4,
     >                      NSR,STRA)
      IF(NSR.GT.I4) CALL ENDJOB('Sbr rrebel. Too large value nsr=',NSR)
      READ(STRA(1),*,ERR=98) A(NOEL,1)
      NRBLT = NINT(A(NOEL,1))
      READ(STRA(2),*,ERR=98) A(NOEL,2)
      READ(STRA(3),*,ERR=98) A(NOEL,3)
      IA3 = NINT(A(NOEL,3))

      IOP = 0
      IF(STRCON(STRA(3),'.',
     >                      II)) THEN
        READ(STRA(3)(II+1:FINSTR(STRA(3))),*,ERR=98,END=98) IOP
      ENDIF

      IF(IOP .NE. 1 .AND. IOP .NE. 2) THEN
        IF(NSR .EQ. I4) THEN
          READ(STRA(4),*,ERR=77) IA4
          GOTO 78
 77       IA4 = 0
 78       CONTINUE
        ENDIF
      ELSE
        IA4 = 0
      ENDIF
      A(NOEL,4) = IA4

      IF    (IA4 .EQ. 1) THEN
C Will 'REBELOTE' using new value (as changed by 'REBELOTE' itself) for parameter #KPRM in element #KLM.

        IF(NRBLT .GT. MXLST) THEN
          WRITE(NRES,*)
     >    'SBR RREBEL. Parameter list too large. Has to be .le. NRBLT= '
     >    ,NRBLT
          GOTO 98
        ENDIF
        LINE = LINE + 1
        READ(NDAT,*,ERR=98,END=98) NPRM
        IF(NPRM .GT. MXPRM) THEN
          WRITE(*,*)
     >    'SBR RREBEL. Too many parameters, has to be .le. ',MXPRM
          GOTO 98
        ENDIF
        A(NOEL,10) = NPRM

        DO IPRM = 1, NPRM
          LINE = LINE + 1
          READ(NDAT,FMT='(A)',ERR=98,END=98) TXT300
          IF(STRCON(TXT300,'    !',
     >                          II))
     >    TXT300 = TXT300(DEBSTR(TXT300):FINSTR(TXT300))

          READ(TXT300,*,ERR=98,END=98) STRING

C Two ways to define the element with parameter to be changed : either its keyword, or its number in the sequence
          IF(ISNUM(STRING)) THEN

            READ(TXT300,*,ERR=79,END=79)
     >        KLM, KPRM, TXT300
            TPRM(IPRM,1) = ' '
            TPRM(IPRM,2) = ' '
            TPRM(IPRM,3) = ' '

          ELSE

            NBTP = 1
            IF(STRCON(TXT300,'[',
     >                            II1)) THEN
              NBTP = NBTP + 1
              IF(STRCON(TXT300,']',
     >                              II2)) THEN
                TPRM(IPRM,1) = TXT300(1:II1-1)
                TPRM(IPRM,1) = TPRM(IPRM,1)
     >          (DEBSTR(TPRM(IPRM,1)):FINSTR(TPRM(IPRM,1)))
                TPRM(IPRM,2) = TXT300(II1+1:II2-1)
                TPRM(IPRM,2) = TPRM(IPRM,2)
     >          (DEBSTR(TPRM(IPRM,2)):FINSTR(TPRM(IPRM,2)))
              ELSE
                WRITE(abs(nres),*)
     >          'Pgm rrebel. Stopped while reading parameter data.'
                WRITE(6,*)
     >          'Pgm rrebel. Stopped while reading parameter data.'
                CALL ENDJOB('Input error in keyword[label] data at'
     >          //' LINE=',LINE)
              ENDIF
            ELSE
              TPRM(IPRM,1) = STRING(DEBSTR(STRING):FINSTR(STRING))
              TPRM(IPRM,2) = ' '
            ENDIF
            TPRM(IPRM,3) = ' '

            IEL = 1
            OKKLE = .FALSE.
            DO WHILE(.NOT. OKKLE .AND. IEL .LE. NOEL)
              IF( TPRM(IPRM,1) .EQ.
     >          KLE(IQ(IEL))(DEBSTR(KLE(IQ(IEL))):FINSTR(KLE(IQ(IEL))))
     >          .AND.
     >          (TPRM(IPRM,2) .EQ. LABEL(IEL,1)
     >           .OR.  TPRM(IPRM,2) .EQ. ' ' ))  OKKLE = .TRUE.
c                 write(*,*) iel,KLE(IQ(IEL)),label(iel,1)
c                 write(*,*) iel,TPRM(IPRM,1) ,TPRM(IPRM,2)
              IEL = IEL + 1
            ENDDO

            IF(OKKLE) THEN
              KLM = IEL - 1
C Keyword with parameter to be changed
C            TPRM(IPRM,1) = STRING(DEBSTR(STRING):FINSTR(STRING))
              BACKSPACE(NDAT)
              LINE = LINE + 1
              READ(NDAT,FMT='(A)',ERR=98,END=98) TXT300
              CALL STRGET(TXT300,I3,
     >                              NSR,STRA)
              TPRM(IPRM,3) = ' '
            ELSE
                WRITE(abs(nres),fmt='(a)')
     >          'Pgm rrebel. Could not find keyword[label] '
     >          // TPRM(IPRM,1) //'['//TPRM(IPRM,2)//']'
     >          //'  in zgoubi data sequence.'
                WRITE(6,fmt='(a)')
     >          'Pgm rrebel. Could not find keyword[label] '
     >          // TPRM(IPRM,1) //'['//TPRM(IPRM,2)//']'
     >          //'  in zgoubi data sequence.'
                 CALL ENDJOB('Check keyword[label] data.',-99)
            ENDIF

          ENDIF

          IF(STRCON(STRA(I3),':',
     >                           IS)) THEN
C DO loop style, V1:V_NRBLT
            READ(TXT300,*,ERR=98,END=98) STRING, KPRM
            READ(STRA(I3)(1:IS-1),*,ERR=98,END=98) VA
            READ(STRA(I3)(IS+1:FINSTR(STRA(I3))),*,ERR=98,END=98) VB
            IF    (NRBLT.LE.0) THEN
            ELSEIF(NRBLT.EQ.1) THEN
              PARAM(IPRM,1) = VA
            ELSEIF(NRBLT.LE.MXLST) THEN
              DV = (VB-VA)/DBLE(NRBLT-1)
              DO I = 1, NRBLT
                PARAM(IPRM,I) = VA + DBLE(I-1)*DV
              ENDDO
            ELSE
              WRITE(NRES,*)
     >        'Sbr rrebel. List too large, must contain '
     >        //'number of data .le. nrblt= ',NRBLT
              GOTO 98
            ENDIF
          ELSE
C List style, V1, V2, ..., V_NRBLT
            BACKSPACE(NDAT)
C            READ(TXT300,*,ERR=98,END=98)
            READ(NDAT,*,ERR=98,END=98)
     >      STRING, KPRM, (PARAM(IPRM,I),I=1,NRBLT)
          ENDIF
          GOTO 80

 79       CONTINUE
          WRITE(6,FMT='(A)') ' '
          WRITE(6,FMT='(A)')
     >    'Sbr rrebel. Stopped while reading list of parameter values.'
     >    //' Give NRBLT .le. number of values in list.'
          GOTO 98

 80       CONTINUE
          A(NOEL,20+10*(IPRM-1)) = KLM
          A(NOEL,21+10*(IPRM-1)) = KPRM
          IF(NRBLT .GT. MXLST) THEN
            WRITE(NRES,*) 'Sbr rrebel. Too many data, # must be .LE.'
     >      //' nrblt= ',NRBLT
            GOTO 98
          ENDIF

        ENDDO

        CALL REBEL4(PARAM,TPRM)
        NOELA = 1
        NOELB = NOEL

      ENDIF

C      IF(STRCON(STRA(3),'.',
C     >                      II)) THEN
C        READ(STRA(3)(II+1:FINSTR(STRA(3))),*,ERR=98,END=98) IOP
        IF(IOP .GE. 1) THEN
          IF(IOP .EQ. 1) THEN
C GET LABEL, DEDUCE RELATED NOEL : MULTI-PASS TRACKING WILL LOOP OVER NOEL-REBELOTE
            READ(TXT300,*,ERR=98,END=98) DUM,DUM,DUM,TXTA
            DO JJ = 1, NOEL
              IF(LABEL(JJ,1).EQ.TXTA) THEN
                NOELA = JJ
                GOTO 12
              ENDIF
            ENDDO
            NOELA = 1
 12         CONTINUE
            NOELB = NOEL
          ELSEIF(IOP .EQ. 2) THEN
C GET 2 LABELS, DEDUCE RELATED NOELS : MULTI-PASS TRACKING WILL LOOP OVER NOEL1-NOEL2
            READ(TXT300,*,ERR=98,END=98) DUM,DUM,DUM,TXTA,TXTB
            DO JJ = 1, NOEL
              IF(LABEL(JJ,1).EQ.TXTA) THEN
                NOELA = JJ
                GOTO 11
              ENDIF
            ENDDO
            NOELA = 1
 11         CONTINUE
            DO JJ = 1, NOEL
              IF(LABEL(JJ,1).EQ.TXTB) THEN
                NOELB = JJ
                GOTO 10
              ENDIF
            ENDDO
            NOELB = NOEL
 10         CONTINUE
          ELSE
            CALL ENDJOB(' Sbr rrebel, No such option IOP=',IOP)
          ENDIF
        ELSE
          NOELA = 1
          NOELB = NOEL
        ENDIF
C      ELSE
C        NOELA = 1
C        NOELB = NOEL
C      ENDIF

C      CALL REBEL6(NOELA, NOELB)
      CALL REBLT6(NOELA, NOELB)

      RETURN

 98   CONTINUE
      CALL ENDJOB('*** Pgm robjet, keyword ''REBELOTE'' : '//
     >'input data error, at line #',LINE)
      RETURN

      END
