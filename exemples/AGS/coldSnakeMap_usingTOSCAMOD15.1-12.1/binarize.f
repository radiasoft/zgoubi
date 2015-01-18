      implicit double precision (a-h,o-z)
      CHARACTER*220 TITL
      CHARACTER*220 fname
      integer debstr, finstr

      write(*,*) ' Give map file name : '
      read(*,*) fname

      write(*,*) ' Give ixma, jyma, kzma : '
      read(*,*) ixma, jyma, kzma

      write(*,*) ' Map file name : '
     >  ,fname(debstr(fname):finstr(fname))

      write(*,*) ' ixma, jyma, kzma : ', ixma, jyma, kzma

       read(*,*)

      lr = 1
      lw = 2
      open(unit=lr,file=fname(debstr(fname):finstr(fname)))
      open(unit=lw,file='binarize.out',FORM='UNFORMATTED')

      nhd = 9   ! 9 line header
              DO II=1, NHD
                READ(lr,FMT='(A120)') TITL
                WRITE(*,*) TITL
                WRITE(lw) TITL
              ENDDO

           DO J=1,JYMA        
             JTC = J
             DO  K = 1,KZMA      
               kzc = k
               DO I=1,IXMA            

                   READ(lr,*,end=99) YH,ZH,XH, 
     >                                  BREAD2,BREAD3,BREAD1
                   write(*,*)YH,ZH,XH,BREAD2,BREAD3,BREAD1
                   write(lw)YH,ZH,XH,BREAD2,BREAD3,BREAD1

               ENDDO
             ENDDO
           ENDDO

 99   continue
      close(lr)
      close(lw)

      open(unit=lw,file='binarize.check.out')
      open(unit=lr,file='binarize.out',FORM='UNFORMATTED')

      nhd = 4   ! 4 line header

              DO II=1, NHD
                READ(lr) TITL
                WRITE(lw,*) TITL
              ENDDO

           DO J=1,JYMA        
             JTC = J
             DO  K = 1,KZMA      
               kzc = k
               DO I=1,IXMA            

                   READ(lr,end=98) YH,ZH,XH, 
     >                                  BREAD2,BREAD3,BREAD1
                   write(lw,fmt='(1p,6(e21.12))')
     >                  YH,ZH,XH,BREAD2,BREAD3,BREAD1

               ENDDO
             ENDDO
           ENDDO

 98   continue
      close(lr)
      close(lw)



      stop
      end

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








