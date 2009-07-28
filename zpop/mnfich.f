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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE MNFICH(TYP) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER * (*) TYP

      INTEGER DEBSTR,FINSTR
      CHARACTER * 80  NOMFIC
      CHARACTER * 80  CMND 
C      CHARACTER REP
                                                               
      GOTO 21

20    CONTINUE      
C      WRITE(6,FMT='(/,''  Press RETURN for more'')') 
C      READ(5,200,ERR=20) REP 
C200   FORMAT(A1)
      I=IDLG('('' Press RETURN for more :'')','    ',1)
 
21    CONTINUE
      CALL HOMCLR
 
22    CONTINUE
 
      WRITE(6,104) (TYP,I=1,6) 
104   FORMAT(//,3X,60(1H*),//,10X,' MENU  - FILES  ',A4,'  :',//
     1 ,9X,'  1        SAVE   (FROM  ZGOUBI', A4 ,'  TO  FNAME',/
     2 ,9X,'  2        READ   (FROM  FNAME  TO  ZGOUBI',A4, /
     3 ,9X,'  3        DELETE ( FNAME ) ',/
     5 ,9X,'  5        DIRECTORY             ',/
     6 ,9X,'  6        EDIT   ZGOUBI',A4,'   ',/
     7 ,9X,'  7        PRINT  ZGOUBI',A4,'   ',/
     8 ,9X,'  8        VISU   ZGOUBI',A4,'   ',/
     9 ,9X,'  9        EXIT THIS MENU        ',///,3X,60(1H*),/)

      WRITE(6,105) ' * Option  number : '
 105  FORMAT(A20,$)
      READ(5,108,ERR=21) IOPT
108   FORMAT(I2)
      GOTO ( 1, 2, 3,20, 5, 6, 7, 8,99) IOPT  
      GOTO 21
 
 1    CONTINUE
         WRITE(6,101)
 101     FORMAT('$ NOM DU FICHIER DE SAUVEGARDE   (QUIT = QQ) : ')
         READ(5,106,ERR=1) NOMFIC
         N1 = DEBSTR(NOMFIC)                            
         N2 = FINSTR(NOMFIC) 
         IF(NOMFIC(N1:N2).EQ. 'QQ' .OR. NOMFIC(N1:N2).EQ. 'qq') GOTO 20
C         CMND= 'COPY zgoubi'//TYP//' '//NOMFIC(N1:N2)//TYP 
         CMND= 'COPY zgoubi'//TYP//' '//NOMFIC(N1:N2)
         WRITE(6,*) '   ',CMND,'EN COURS...' 
         IOUPI=LIBSPA(CMND)  
         WRITE(6,*) '   ',CMND,'TERMINEE...' 
      GOTO 20                

 2    CONTINUE
         CALL DRCTRY(TYP)
         WRITE(6,102)
 102     FORMAT('$ NOM DU FICHIER A LIRE (QUIT = QQ) : ')
         READ(5,106,ERR=2) NOMFIC
106      FORMAT(A50)
         N1 = DEBSTR(NOMFIC)                            
         N2 = FINSTR(NOMFIC) 
         IF(NOMFIC(N1:N2).EQ. 'QQ' .OR. NOMFIC(N1:N2).EQ. 'qq') GOTO 20
C         CMND= 'COPY '//NOMFIC(N1:N2)//TYP//' zgoubi'//TYP  
         CMND= 'COPY '//NOMFIC(N1:N2)//' zgoubi'//TYP  
         WRITE(6,*) '   ',CMND,'EN COURS...' 
         IOUPI=LIBSPA(CMND)
         WRITE(6,*) '   ',CMND,'TERMINEE...' 
      GOTO 20                

 3    CONTINUE
         CALL DRCTRY(TYP)
         WRITE(6,103)
 103     FORMAT('$ NOM DU FICHIER A DELETER (QUIT = QQ) : ')
         READ(5,106,ERR=3) NOMFIC       
         N1 = DEBSTR(NOMFIC)                            
         N2 = FINSTR(NOMFIC) 
         IF(NOMFIC(N1:N2).EQ. 'QQ' .OR. NOMFIC(N1:N2).EQ. 'qq') GOTO 20
         CMND= 'DELETE/ERASE/CONFIRM '//NOMFIC(N1:N2)
         WRITE(6,*) '   ',CMND,'EN COURS...' 
         IOUPI=LIBSPA(CMND)
         WRITE(6,*) '   ',CMND,'TERMINEE...' 
      GOTO 20                

C         CALL DRCTRY(TYP)
C         WRITE(6,105)
C 105     FORMAT('$ NOM DU FICHIER A RENAMER (QUIT = QQ) : ')
C         READ(5,106,ERR=4) NOMFIC                       
C         N1 = DEBSTR(NOMFIC)                            
C         N2 = FINSTR(NOMFIC) 
C         IF(NOMFIC(N1:N2).EQ. 'QQ' .OR. NOMFIC(N1:N2).EQ. 'qq') GOTO 20
C         CMND= 'RENAME '//NOMFIC(N1:N2)//TYP 
C         IOUPI=LIBSPA(CMND)
C         WRITE(6,*) '   ',CMND,'TERMINEE...' 
      GOTO 20 
 5    CONTINUE                  
         CALL DRCTRY(TYP)
      GOTO 20 
 6    CONTINUE                                          
         CALL HOMCLR 
         IF(TYP .NE. '.dat') IOUPI=LIBSPA('SET TERM/WIDTH=132') 
         CMND= 'ED zgoubi'//TYP 
         WRITE(6,*) '   ',CMND,'EN COURS...' 
         IOUPI=LIBSPA(CMND) 
         CALL HOMCLR 
         IF(TYP .NE. '.dat') IOUPI=LIBSPA('SET TERM/WIDTH=80 ') 
      GOTO 22
 7    CONTINUE  
C         CMND='PR/HEAD/CONF/QUE=ERICH$LASER_HH zgoubi'//TYP//' '
         CMND=' HPRINT  zgoubi'//TYP//' E '
         IOUPI=LIBSPA(CMND) 
      GOTO 20       
 8    CONTINUE
         CALL HOMCLR 
         IF(TYP .NE. '.dat') IOUPI=LIBSPA('SET TERM/WIDTH=132') 
         CMND='COPY zgoubi'//TYP//' SYS$OUTPUT'  
C         CMND='MORE zgoubi'//TYP
         WRITE(6,*) '   ',CMND,'EN COURS...' 
         IOUPI=LIBSPA(CMND)  
C         WRITE(6,FMT='(/,''  Press RETURN for more'')') 
C         READ(5,200) REP 
         I=IDLG('('' Press RETURN for more :'')','    ',1)
         CALL HOMCLR 
         IF(TYP .NE. '.dat') IOUPI=LIBSPA('SET TERM/WIDTH=80 ') 
      GOTO 22       
C                                  
99    RETURN
      END 
