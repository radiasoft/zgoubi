c     positi.f
c
c     written by Frédéric Desforges, frederic.desforges@grenoble-inp.org
c
c
      SUBROUTINE POSITI(NRES,ARCLEN,LABEL,KEYWOR)

      IMPLICIT DOUBLE PRECISION (A-H,M-Z)
      INTEGER NRES,IS,I,DEBSTR,FINSTR
      CHARACTER*300 BUFFER,LABEL,KEYWOR
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
