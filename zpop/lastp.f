      SUBROUTINE LASTP(NL,
     >                    NPASS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C---------------------------------
C Look for last pass number, NPASS
C---------------------------------
      CHARACTER LET 

      INCLUDE 'MXVAR.H'
      DIMENSION YZXB(MXVAR),NDX(5)

          CALL REWIN2(NL,*96)
          WRITE(6,*) '  Reading coordinates, looking for NPASS, wait...'
          NOC=0
          NRBLT = -1 
C--------- BOUCLE SUR READ FICHIER NL 
 45       CONTINUE
            CALL READCO(NL,
     >                        KART,LET,YZXB,NDX,*12,*78)

C--------- NDX: 1->KEX, 2->IT, 3->IREP, 4->IMAX

            NOC=NOC+1

            IF(NINT(YZXB(39)) .GE. NRBLT+1) NRBLT = NINT(YZXB(39)) -1
            IF(NOC.EQ. NPTR) GOTO 13

          GOTO 45
C         ----------------------------------
 78       CONTINUE
          WRITE(6,*) ' *** Coordinates reading stopped : error during',
     >     ' read of event # ',NOC+1
          GOTO 13

 12       CONTINUE
          WRITE(6,*) ' READ  OK; END  OF  FILE  ENCOUNTERED'

 13       CONTINUE
          NPASS = NRBLT + 1
          WRITE(6,*) ' Found last pass number, NPASS = ',NPASS

 96   RETURN
      END
