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
      CHARACTER(I300) TXT300, TX300
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
      LOGICAL ISNUM, EMPTY
      PARAMETER (I3=3)

      CHARACTER(10) TXT10

      DATA IA4 / 0 /

      LINE = 1
      READ(NDAT,FMT='(A)',ERR=98,END=98) TXT300
      TXT300 = TXT300(DEBSTR(TXT300):FINSTR(TXT300))
      IF(STRCON(TXT300,'        !',
     >                             II))
     >TXT300 = TXT300(DEBSTR(TXT300):II-1)
      CALL STRGET(TXt300,I4,
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
        READ(NDAT,*,ERR=98,END=98) A(NOEL,10)
        NPRM = NINT(A(NOEL,10))
        IF(NPRM .GT. MXPRM) THEN
          WRITE(6,*)
     >    'SBR RREBEL. Too many parameters, has to be .le. ',MXPRM
          GOTO 98
        ENDIF

        DO IPRM = 1, NPRM
          LINE = LINE + 1
          READ(NDAT,FMT='(A)',ERR=98,END=98) TXT300
          IF(STRCON(TXT300,'    !',
     >                             II))
     >    TXT300 = TXT300(DEBSTR(TXT300):FINSTR(TXT300))
          READ(TXT300,*,ERR=98,END=98) STRING
C Two ways to define the element with parameter to be changed : either its keyword, or its number in the sequence
          IF(ISNUM(STRING)) THEN

            READ(TXT300,*,ERR=79,END=79) KLM, KPRM
            TPRM(IPRM,1) = ' '  ! Keyword is absent, element concerned is specified by its number in the sequence
            TPRM(IPRM,2) = ' '  ! label1 absent
            TPRM(IPRM,3) = ' '  ! label2 absent

            WRITE(TXT10,FMT='(I0)') KLM
            TXT300 = TXT300(LEN(TRIM(TXT10))+1:)
            II2 = 0
c            write(*,*) ' ||| txt300 ',trim(txt300)
c            write(*,*) ' *** txt10 ',trim(txt10),II2

          ELSE

            IF(STRCON(TXT300,'{',
     >                            II1)) THEN
              IF(STRCON(TXT300,'}',
     >                              II2)) THEN
                TPRM(IPRM,1) = TXT300(1:II1-1)                             ! Keyword
                TPRM(IPRM,1) = TPRM(IPRM,1)(DEBSTR(TPRM(IPRM,1))
     :          :FINSTR(TPRM(IPRM,1)))

                STRING =TXT300(II1+1:II2-1)

                IF(.NOT. EMPTY(STRING)) THEN
                  IF(STRCON(TXT300,',',
     >                                 II3)) THEN
                    IF(.NOT. EMPTY(TXT300(II1+1:II3-1)) ) THEN
                      READ(TXT300(II1+1:II3-1),*,ERR=79,END=79)
     >                TPRM(IPRM,2)                                          ! LABEL1
                      TPRM(IPRM,2) = TPRM(IPRM,2)
     >                (DEBSTR(TPRM(IPRM,2)):FINSTR(TPRM(IPRM,2)))
                    ELSE
                      TPRM(IPRM,2) = ' '
                    ENDIF
                    IF(.NOT. EMPTY(TXT300(II3+1:II2-1))) THEN
                      READ( TXT300(II3+1:II2-1),*,ERR=79,END=79)
     >                TPRM(IPRM,3)                                         ! LABEL2
                      TPRM(IPRM,3) = TPRM(IPRM,3)
     >                (DEBSTR(TPRM(IPRM,3)):FINSTR(TPRM(IPRM,3)))
                    ELSE
                      TPRM(IPRM,3) = ' '
                    ENDIF
                  ELSE
                    READ( TXT300(II1+1:II2-1),*,ERR=79,END=79)
     >              TPRM(IPRM,2)                                          ! LABEL1
                    TPRM(IPRM,2) = TPRM(IPRM,2)
     >              (DEBSTR(TPRM(IPRM,2)):FINSTR(TPRM(IPRM,2)))
                    TPRM(IPRM,3) = ' '
                  ENDIF

                ELSE
                  TPRM(IPRM,2) = ' '
                  TPRM(IPRM,3) = ' '
                ENDIF

              ELSE
                WRITE(ABS(NRES),*)
     >          'Pgm rrebel. Stopped while reading parameter data.'
     >          //' Missing ''}'' ?'
                WRITE(6,*)
     >          'Pgm rrebel. Stopped while reading parameter data.'
     >          //' Missing ''}'' ?'
                CALL ENDJOB('Input error in keyword{label} data at'
     >          //' LINE=',LINE)
              ENDIF
            ELSE
              TPRM(IPRM,1) = STRING(DEBSTR(STRING):)
              TPRM(IPRM,2) = ' '   ! No LABEL1
              TPRM(IPRM,3) = ' '   ! No LABEL2
              II2 = LEN(TRIM(STRING))
            ENDIF

            IEL = 1
            OKKLE = .FALSE.
            DO WHILE(.NOT. OKKLE .AND. IEL .LE. NOEL)
              IF( (TRIM(TPRM(IPRM,1)) .EQ. TRIM(KLE(IQ(IEL))))
     >        .AND.
     >        ((TRIM(TPRM(IPRM,2)) .EQ. TRIM(LABEL(IEL,1)) )
     >        .OR.  EMPTY(TPRM(IPRM,2)))
     >        .AND.
     >        ((TRIM(TPRM(IPRM,3)) .EQ. TRIM(LABEL(IEL,2)) )
     >        .OR.  EMPTY(TPRM(IPRM,3)))
     >         )  OKKLE = .TRUE.
              IEL = IEL + 1
            ENDDO

            IF(OKKLE) THEN
              KLM = IEL - 1
C Keyword with parameter to be changed
              BACKSPACE(NDAT)
              LINE = LINE + 1
              READ(NDAT,FMT='(A)',ERR=98,END=98) TXT300
              TX300 = TXT300
              CALL STRGET(TX300,I3,
     >                             NSR,STRA)
C              TPRM(IPRM,3) = ' '

            ELSE
                WRITE(abs(nres),fmt='(a)')
     >          'Pgm rrebel. Could not find keyword{label} '
     >          // TPRM(IPRM,1) //'{'//TPRM(IPRM,2)//'}'
     >          //'  in zgoubi data sequence.'
                WRITE(6,fmt='(a)')
     >          'Pgm rrebel. Could not find keyword{label} '
     >          // TPRM(IPRM,1) //'{'//TPRM(IPRM,2)//'}'
     >          //'  in zgoubi data sequence.'
                 CALL ENDJOB('Check keyword{label} data.',-99)
            ENDIF

          ENDIF

c             write(*,*) ' ii2 = ',iI2
c             write(*,*) ' string ',iprm,trim(string),' ',TPRM(IPRM,1)
c     >              ,' * ',TPRM(IPRM,2),' * ',TPRM(IPRM,3),' ',okkle
c             write(*,*) ' TXT300 ',TXT300
c              write(*,*) ' TXT300(II2+1:) ',trim(TXT300(II2+1:))
c              read(*,*)

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
C            BACKSPACE(NDAT)

C            IF(STRCON(STRING,',',
C     >                           IS4)) THEN
C     READ(NDAT,*,ERR=98,END=98)

C                write(*,*) ' length(string)= ',len(trim(string))

              READ(TXT300(II2+1:),*,ERR=98,END=98)
C              READ(TXT300(len(trim(string))+1:),*,ERR=98,END=98)
     >        KPRM, (PARAM(IPRM,I),I=1,NRBLT)
C            ELSE
C              READ(NDAT,*,ERR=98,END=98)
C     >        STRING, KPRM, (PARAM(IPRM,I),I=1,NRBLT)
C            ENDIF
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

      CALL REBLT6(NOELA, NOELB)

      RETURN

 98   CONTINUE
      CALL ENDJOB('*** Pgm rrebel, keyword ''REBELOTE'' : '//
     >'input data error, at line #',LINE)
      RETURN

      END
