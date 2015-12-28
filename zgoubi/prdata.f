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
      SUBROUTINE PRDATA(NLIN,FLIN,FDAT,NRES,
     >                                 LABEL,NOEL,NDAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) FLIN, FDAT
      INCLUDE 'MXLD.H'
      CHARACTER(*) LABEL(MXL,*)
C----------------------------------------------------------
C     Copy zgoubi.dat into zgoubi.res, and a few more stuff
C----------------------------------------------------------
      PARAMETER (I2000=2000)
      CHARACTER(I2000) TEXT
      PARAMETER (I6=6)
      CHARACTER(I6) TXT6
      INTEGER DEBSTR, FINSTR
      LOGICAL EMPTY
      LOGICAL IDLUNI, OK

      PARAMETER (LBLSIZ=10)
      CHARACTER(LBLSIZ) LAB2(2)
      PARAMETER (KSIZ=10)
      PARAMETER (I104=104)
      CHARACTER(110) TXT110
      DIMENSION LUNR(10)
      PARAMETER(MXFIL=1)
      CHARACTER(132) FINC(MXFIL)
      CHARACTER(50) CMMND
      LOGICAL YINC, YINC2, STRCON
      CHARACTER(LBLSIZ) LBL1A(MXFIL),LBL2A(MXFIL),
     >LBL1B(MXFIL),LBL2B(MXFIL)
      CHARACTER(LBLSIZ) L1A,L2A,L1B,L2B 
      LOGICAL LBAVU, LBBVU, EXS

      DATA YINC, YINC2 / .FALSE. , .FALSE. /
      DATA EXS / .FALSE. /

      WRITE(6,*) '  Copying  zgoubi.dat  into  zgoubi.res,'
      WRITE(6,*) '  numbering  and  labeling  elements...'

      IDA = 1
      LUNR(IDA) = NLIN
      YINC = .FALSE. 
      YINC2 = .FALSE. 
      
C----- Read zgoubi.dat title (1st data line)
      READ(LUNR(IDA),FMT='(A)',ERR=10,END=95) TEXT
      WRITE(NRES,FMT='(A)') TEXT(DEBSTR(TEXT):FINSTR(TEXT))

      NOEL=0
      L1A = '*'      
      L2A = '*'      
      L1B = '*'      
      L2B = '*'      
 10   CONTINUE
        READ (LUNR(IDA),FMT='(A)',ERR=10,END=95) TEXT
        IF( .NOT. EMPTY(TEXT) ) THEN
          TEXT = TEXT(DEBSTR(TEXT):FINSTR(TEXT))
        ELSE
          TEXT = ' '
        ENDIF

        IDEB = 1
        IF( TEXT(IDEB:IDEB) .EQ. '''' ) THEN
          NOEL=NOEL+1
          TEXT = TEXT(DEBSTR(TEXT):I104)

          DO I=IDEB+1,IDEB+KSIZ+1
            IF(TEXT(I:I) .EQ. '''') GOTO 2
          ENDDO

 2        CONTINUE
          
          LABEL(NOEL,1) = ' '
          LABEL(NOEL,2) = ' '
          IF( .NOT. EMPTY(TEXT((I+1):I2000)) ) THEN
            CALL STRGET(TEXT((I+1):I2000),2,
     >                                      NST,LAB2)
            IF(NST.GE.1) THEN
              IF(LAB2(1)(1:1).NE.'!') THEN
                IF(.NOT. EMPTY(LAB2(1))) LABEL(NOEL,1) = LAB2(1)
                IF(NST.EQ.2) THEN
                  IF(LAB2(2)(1:1).NE.'!') THEN
                    IF(.NOT. EMPTY(LAB2(2))) LABEL(NOEL,2) = LAB2(2)
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF

          IF(YINC2) THEN
            IF(LBAVU) THEN
              LBBVU = ((L1B .EQ. '*') .AND. (L2B .EQ. LABEL(NOEL,2)))
     >        .OR.    ((L1B .EQ. LABEL(NOEL,1)) .AND. (L2B .EQ. '*'))
     >        .OR.((L1B.EQ.LABEL(NOEL,1)).AND.(L2B.EQ.LABEL(NOEL,2)))

              IF(LBBVU) THEN
                WRITE(TEXT,FMT='(A)') TEXT(DEBSTR(TEXT):FINSTR(TEXT))
     >          //'  ! Include_End : '
     >          //FINC(IFL)(DEBSTR(FINC(IFL)):FINSTR(FINC(IFL)))
                WRITE(TXT6,FMT='(I6)') NOEL
                TEXT = TEXT(1:I104)//TXT6
                WRITE(NRES,FMT='(T2,A)') TEXT(1:110)
              ENDIF

            ELSE
              LBAVU = ((L1A .EQ. '*') .AND.
     >        ((L2A .EQ. '*') .OR. (L2A .EQ. LABEL(NOEL,2))))
     >        .OR.    ((L1A .EQ. LABEL(NOEL,1)) .AND.
     >        ((L2A .EQ. '*') .OR. (L2A .EQ. LABEL(NOEL,2))))

              IF(LBAVU) THEN
                IF((L1A .EQ. '*') .AND. (L2A .EQ. '*')) THEN 
                  WRITE(TXT110,FMT='(A)') '''MARKER'''
     >            //'  ! Include_Start :'
                ELSE
                  WRITE(TXT110,FMT='(A)')TEXT(DEBSTR(TEXT):FINSTR(TEXT))
     >            //'  ! Include_Start :'
                ENDIF
                WRITE(TXT110,FMT='(A)') 
     >          TXT110(DEBSTR(TXT110):FINSTR(TXT110))
     >          //' '//FINC(IFL)(DEBSTR(FINC(IFL)):FINSTR(FINC(IFL)))
     >          //';  range : '
     >          //'['//LBL1A(IFL)(DEBSTR(LBL1A(IFL)):FINSTR(LBL1A(IFL)))
     >          //','//LBL2A(IFL)(DEBSTR(LBL2A(IFL)):FINSTR(LBL2A(IFL)))
     >          //':'//LBL1B(IFL)(DEBSTR(LBL1B(IFL)):FINSTR(LBL1B(IFL)))
     >          //','//LBL2B(IFL)(DEBSTR(LBL2B(IFL)):FINSTR(LBL2B(IFL)))
     >          //']'
                WRITE(TXT6,FMT='(I6)') NOEL
                TXT110 = TXT110(1:I104)//TXT6
                WRITE(NRES,FMT='(T2,A)') TXT110 
                if((l1a .eq. '*') .and. (l2a .eq. '*')) then 
                  NOEL = NOEL + 1
                  WRITE(TXT6,FMT='(I6)') NOEL
                  TXT110 = TEXT(1:I104)//TXT6
                  WRITE(NRES,FMT='(T2,A)')
     >            TXT110(DEBSTR(TXT110):FINSTR(TXT110))
                ENDIF
                GOTO 10
              ENDIF
            ENDIF
            IF( .NOT. LBAVU ) THEN
              NOEL = NOEL -1 
              GOTO 10
            ENDIF
            IF(  LBBVU ) GOTO 95
          ENDIF

          IF(TEXT(IDEB+1:I-1) .EQ. 'INCLUDE') THEN
             
            YINC2 = .TRUE.

            IF(FLIN(FINSTR(FLIN)-9:FINSTR(FLIN)) .EQ. 'zgoubi.dat')
     >      CALL ENDJOB('Pgm prdata. Job includes INCLUDE keyword, '//
     >      'hence "zgoubi.dat" cannot be used as input file name '//
     >      '(zgoubi.dat is reserved). Please use different input '
     >      //'file name. Use "zgoubi -fileIn filename" command.'
     >      ,-99)

            LINE =1 
            READ(LUNR(IDA),FMT='(A)',ERR=10,END=95) TEXT
            READ(TEXT,*,ERR=78) NBFIL
            IF(NBFIL .GT. 1) CALL ENDJOB('Pgm prdata. '//
     >      'INCLUDE is limited to 1 file, '//
     >      '(functionality of NBFIL>1 needs be checked). '//
     >      'Use several INCLUDEs instead.',-99)

            IF(NBFIL.GT.MXFIL) CALL ENDJOB('Pgm prdata. '//
     >      'INCLUDE has too many files. Max allowed is ',MXFIL)
            DO IFL = 1, NBFIL
              LINE = LINE + 1 
              READ(LUNR(IDA),FMT='(A)',ERR=10,END=95) TEXT
C TEXT is of the form FILENAME[LBL1a,LBLl2a:LBL1b,LBL2b]. 
C Only text within range Any_KEYWORD LABEL1=lbl1a&LABEL2=lbl2a : Any_KEYWORD LABEL1=lbl1b&LABEL2=lbl2b 
C will be included
              OK = STRCON(TEXT,'[',
     >                             ISA) 
              IF(OK) THEN
                FINC(IFL) = TEXT(DEBSTR(TEXT):ISA-1)
              ELSE
                FINC(IFL) = TEXT(DEBSTR(TEXT):FINSTR(TEXT))
              ENDIF

              IF(OK) THEN
                OK = STRCON(TEXT,']',
     >                                     ISB) 
                IF(.NOT. OK) CALL ENDJOB
     >          ('Sbr prdata. Keyword INCLUDE. Formatting error '//
     >          ' in INCLUDE[lbl1a,lbl2a:lbl1b,lbl2b] (missing '//
     >          ' closing square braket ?)',-99 )

                OK = STRCON(TEXT,':',
     >                               ISC)
                IF(.NOT. OK) CALL ENDJOB
     >          ('Sbr prdata. Keyword INCLUDE. Formatting error '//
     >          ' in INCLUDE[lbl1a,lbl2a:lbl1b,lbl2b] (missing '//
     >          ' ":" ?)',-99 )

                OK = STRCON(TEXT(1:ISC-1),',',
     >                                        ISAC)

                IF(OK) THEN
                
                  IF(ISA+1 .LE. ISAC-1) THEN
C TEXT is of the form FILENAME[lbl1a,lbll2a: NOT YET KNOWN]
                    LBL1A(IFL) = TEXT(ISA+1:ISAC-1)
                  ELSE
                    LBL1A(IFL) = '*'
                  ENDIF

                  IF(ISAC+1 .LE. ISC-1) THEN
                    LBL2A(IFL) = TEXT(ISAC+1:ISC-1)
                  ELSE
                    LBL2A(IFL) = '*'
                  ENDIF

                ELSE

                  IF(ISA+1 .LE. ISC-1) THEN
                    LBL1A(IFL) = TEXT(ISA+1:ISC-1)
                  ELSE
                    LBL1A(IFL) = '*'
                  ENDIF
                  LBL2A(IFL) = '*'
                ENDIF

                OK = STRCON(TEXT(ISC+1:ISB-1),',',
     >                                          ISCB)
                ISCB = ISCB + ISC

                IF(OK) THEN

                  IF(ISC+1 .LE. ISCB-1) THEN
                    LBL1B(IFL) = TEXT(ISC+1:ISCB-1)
                  ELSE
                    LBL1B(IFL) = '*'
                  ENDIF

                  IF(ISCB+1 .LE. ISB-1) THEN
                    LBL2B(IFL) = TEXT(ISCB+1:ISB-1)
                  ELSE
                    LBL2B(IFL) = '*'
                  ENDIF

                ELSE

                  IF(ISC+1 .LE. ISB-1) THEN
                    LBL1B(IFL) = TEXT(ISC+1:ISB-1)
                  ELSE
                    LBL1B(IFL) = '*'
                  ENDIF
                  LBL2B(IFL) = '*'
                ENDIF

                LBL1A(IFL) = 
     >          LBL1A(IFL)(DEBSTR(LBL1A(IFL)):FINSTR(LBL1A(IFL)))
                LBL2A(IFL) = 
     >          LBL2A(IFL)(DEBSTR(LBL2A(IFL)):FINSTR(LBL2A(IFL)))
                LBL1B(IFL) = 
     >          LBL1B(IFL)(DEBSTR(LBL1B(IFL)):FINSTR(LBL1B(IFL)))
                LBL2B(IFL) = 
     >          LBL2B(IFL)(DEBSTR(LBL2B(IFL)):FINSTR(LBL2B(IFL)))

              ELSE

                LBL1A(IFL) = '*'
                LBL2A(IFL) = '*'
                LBL1B(IFL) = '*'
                LBL2B(IFL) = '*'

              ENDIF

            ENDDO
            IDA = IDA + 1
            OK = IDLUNI(
     >                  LUNR(IDA)) 
            IFL = 1

            INQUIRE(FILE=FINC(IFL),EXIST=EXS)
            IF(.NOT. EXS) CALL ENDJOB('Pgm prdata, keyword INCLUDE  : '
     >      //'could not include file, does not exist ',-99)
            OPEN(UNIT=LUNR(IDA),FILE=FINC(IFL))           

            L1A = LBL1A(IFL)
            L2A = LBL2A(IFL)
            L1B = LBL1B(IFL)
            L2B = LBL2B(IFL)

            NOEL = NOEL -1 
            GOTO 10

          ELSE

            WRITE(TXT6,FMT='(I6)') NOEL
            TEXT = TEXT(1:I104)//TXT6
            WRITE(NRES,FMT='(T2,A)') TEXT(1:110)
          ENDIF

          IF(   TEXT(IDEB:IDEB+4) .EQ. '''FIN'''
     >     .OR. TEXT(IDEB:IDEB+4) .EQ. '''END''') GOTO 95

        ELSE

            if( yinc2 .and. ((.not. lbavu) .or. lbbvu )) goto 10

          WRITE(NRES,FMT='(A)') TEXT(DEBSTR(TEXT):FINSTR(TEXT))
c          WRITE(NTMP,FMT='(A)') TEXT(DEBSTR(TEXT):FINSTR(TEXT))

        ENDIF

      GOTO 10

 95   CONTINUE

      IF(IDA.GT.1) THEN
        IF(YINC2) THEN 
          IF(.NOT. LBBVU) THEN
            WRITE(TXT110,FMT='(A)') '''MARKER'''
     >      //'  ! Include_End : '
     >      //FINC(IFL)(DEBSTR(FINC(IFL)):FINSTR(FINC(IFL)))
            WRITE(TXT6,FMT='(I6)') NOEL
            TXT110 = TXT110(1:I104)//TXT6
            WRITE(NRES,FMT='(T2,A)') TXT110 
          ENDIF
        ENDIF
        YINC2 = .FALSE.
        LBAVU=.FALSE.
        LBBVU=.FALSE.
        YINC = .TRUE.
        CLOSE(LUNR(IDA))
        IDA = IDA-1
        GOTO 10
      ENDIF

      IF(YINC) THEN
        CLOSE(NLIN)
        CALL FLUSH2(NRES,.FALSE.)
        CMMND = 'cp zgoubi.res '//FDAT(DEBSTR(FDAT):FINSTR(FDAT))
        CALL SYSTEM(CMMND)
        OK=IDLUNI(
     >            NDAT)
        OPEN(UNIT=NDAT,FILE=FDAT)
      ELSE
        REWIND(NLIN)
        NDAT = NLIN
      ENDIF

      RETURN

 78   CONTINUE
      CALL ENDJOB('*** Pgm prdata, keyword INCLUDE : '// 
     >'input data error, at line ',LINE)
      RETURN
      END
