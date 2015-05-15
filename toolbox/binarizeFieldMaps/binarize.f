      implicit double precision (a-h,o-z)
C      CHARACTER*120 HDR(88)
      INCLUDE "MAPHDR.H"
      CHARACTER(220) FNAME
      INTEGER DEBSTR, FINSTR
      LOGICAL FIRST

      WRITE(*,*) ' Give map file name : '
      READ(*,*) FNAME

      WRITE(*,*) ' Give ixma, jyma, kzma : '
      READ(*,*) IXMA, JYMA, KZMA

      NHD = -99
      DOWHILE (NHD .LT.0 .OR. NHD .GT.MXHD)
        WRITE(*,*) ' Give number of header lines  (.le. ',
     >  MXHD,')  : '
        READ(*,*) NHD
      ENDDO

      WRITE(*,*) ' Map file name : '
     >  ,FNAME(DEBSTR(FNAME):FINSTR(FNAME))
      WRITE(*,*) ' ixma, jyma, kzma : ', IXMA, JYMA, KZMA
      WRITE(*,*) ' Number of header lines : ', NHD
      WRITE(*,*) ' Ok ? '
      READ(*,*)

      LR = 1
      LW = 2
      LRFN = 3
      OPEN(UNIT=LR,FILE=FNAME(DEBSTR(FNAME):FINSTR(FNAME)))
      OPEN(UNIT=LW,FILE='binarize.out',FORM='UNFORMATTED')

c      nhd = 9   ! 9 line header
           DO II=1, NHD
             READ(lr,FMT='(A)') HDR(II)
             WRITE(*,*) HDR(II)
             WRITE(lw) HDR(II)
           ENDDO

           FIRST = .TRUE.
           DO J=1,JYMA        
             JTC = J
             DO  K = 1,KZMA      
               KZC = K
               DO I=1,IXMA            

                 READ(LR,*,end=99) YH,ZH,XH, 
     >                                  BREAD2,BREAD3,BREAD1
                 IF(FIRST)
     >           WRITE(*,FMT='(A)') 
     >           (HDR(II)(DEBSTR(HDR(II)):FINSTR(HDR(II)))
     >           ,II=1,NHD) 
                                     
C                 write(*,*)YH,ZH,XH,BREAD2,BREAD3,BREAD1
                 WRITE(LW)YH,ZH,XH,BREAD2,BREAD3,BREAD1

                 IF(FIRST) THEN
                   FIRST=.FALSE.
                   WRITE(*,*)YH,ZH,XH,BREAD2,BREAD3,BREAD1
                   WRITE(*,*) 'Ok with ',NHD,'-line header '
     >             //' and first line of field data ?'
                   READ(*,*)

                   WRITE(*,*) ' Busy... working ...'
                 ENDIF

               ENDDO
             ENDDO
           ENDDO

 99   CONTINUE
      CLOSE(LR)
      CLOSE(LW)

      OPEN(UNIT=LW,FILE='binarize.check.out')
      OPEN(UNIT=LR,FILE='binarize.out',FORM='UNFORMATTED')
      OPEN(UNIT=LRFN,FILE=FNAME(DEBSTR(FNAME):FINSTR(FNAME)))

C      nhd = 4   ! 4 line header

              DO II=1, NHD
                READ(LR) HDR(II)
                WRITE(LW,*) HDR(II)(DEBSTR(HDR(II)):FINSTR(HDR(II)))
                READ(LRFN,FMT='(A)') HDR(II)
                WRITE(LW,*) HDR(II)(DEBSTR(HDR(II)):FINSTR(HDR(II)))
              ENDDO

           SUM = 0.D0
           DO J=1,JYMA        
             JTC = J
             DO  K = 1,KZMA      
               KZC = K
               DO I=1,IXMA            

                   READ(LR,END=98) YH,ZH,XH, 
     >                                  BREAD2,BREAD3,BREAD1
                   READ(LRFN,*,END=98) YHi,ZHi,XHi, 
     >                                  BREAD2i,BREAD3i,BREAD1i
                   WRITE(LW,FMT='(1P,6(E21.12))')
     >                  YH,ZH,XH,BREAD2,BREAD3,BREAD1
                   SUM = SUM+ ABS(YH-YHI)+ABS(ZH-ZHI)+ABS(XH-XHI)
     >             + ABS(BREAD1-BREAD1I) + ABS(BREAD2-BREAD2I)  
     >              + ABS(BREAD3-BREAD3I)
               ENDDO
             ENDDO
           ENDDO

         WRITE(*,*) '  '
         WRITE(*,*) '  '
         WRITE(*,*) ' This is the difference between input file and '
     >   //' binary.out. Should be zero. If not then something '
     >   //'went wrong, check !'
         WRITE(*,*) ' difference = ',SUM
         WRITE(*,*) '  '

 98   CONTINUE
      CLOSE(LR)
      CLOSE(LW)



      STOP
      END

      FUNCTION DEBSTR(STRING)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DEBSTR
      CHARACTER * (*) STRING

C     --------------------------------------
C     RENVOIE DANS DEBSTR LE RANG DU
C     1-ER CHARACTER NON BLANC DE STRING,
C     OU BIEN 0 SI STRING EST VIDE ou BLANC.
C     --------------------------------------

      DEBSTR=0
      LENGTH=LEN(STRING)

      IF(LENGTH.EQ.0) RETURN

1     CONTINUE
        DEBSTR=DEBSTR+1
        IF (STRING(DEBSTR:DEBSTR) .EQ. ' ') THEN
          IF(DEBSTR .GE. LENGTH) THEN
            DEBSTR = 0
            RETURN
          ELSE
            GOTO 1
          ENDIF
        ENDIF

      RETURN
      END

      FUNCTION FINSTR(STRING)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER FINSTR
      CHARACTER * (*) STRING
C     --------------------------------------
C     RENVOIE DANS FINSTR LE RANG DU
C     DERNIER CHARACTER NON BLANC DE STRING,
C     OU BIEN 0 SI STRING EST VIDE ou BLANC.
C     --------------------------------------

      FINSTR=LEN(STRING)+1
1     CONTINUE
        FINSTR=FINSTR-1
        IF(FINSTR .EQ. 0) RETURN
        IF (STRING(FINSTR:FINSTR) .EQ. ' ') GOTO 1

      RETURN
      END








