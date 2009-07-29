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
CDECK GRAPH2 
      SUBROUTINE GRAPH8(NLOG, LM, NOMFIC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*)  NOMFIC 
      
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/ OKECH, OKVAR, OKBIN
      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*7, KPOL(2)*9, KDIM(MXVAR)*7
      COMMON/INPVR/ KVAR, KPOL, KDIM
      COMMON/LUN/NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      CHARACTER SAVFIL*20
      LOGICAL OKOPN, INPECH, CHANGE
      SAVE NL
      DATA IZERO /0/
      DATA  OKOPN / .FALSE. /
      DATA CHANGE / .TRUE./
      DATA SAVFIL / 'zpop.ps' /
 
      GOTO 921

 920  CONTINUE      
      CALL FBGTXT
      I=IDLG('('' Press RETURN for more :'')','    ',1)

 921  CONTINUE
      CALL FBGTXT
      CALL HOMCLR

      WRITE(6,100) NOMFIC
 100  FORMAT(//,3X,60(1H*),//,20X,' MENU - Analysis/Graphic :' ,//
     > ,9X,'        ' ,/
     1 ,9X,'  1    Open  file - current is ',A ,/
     3 ,9X,'  3    Plot  options',/
     4 ,9X,'  4    Manual  scales ',/
C     5 ,9X,'  5    Automatic  scales ',/
     6 ,9X,'  6    Synoptic ',/
     8 ,9X,'  8    Save  graphic  screen ',/
     9 ,9X,'  9    Exit  this  menu     ',/
     X ,9X,' 10    SR loss parameters ',/
     1 ,9X,' 11    Isofield  lines from a map',/
     2 ,9X,' 12    Erase  display    ',/
     3 ,9X,' 13    Fourier  spectrum  -  multiturn ',/
     4 ,9X,' 14    Ellipse  fit - multiturn or multiparticle  ',/
     5 ,9X,' 15    Kill/create  graphic  sub-process      ',/
     5 ,9X,' 16    Synchrotron  radiation distributions ',/
     5 ,9X,' 17    Mean  orbit - multiturn c.o., or multiparticle  ',/
     5 ,9X,' 18    Fringe  field  matching',/
     5 ,9X,' 19    Transverse  beam  density  matching',/
     5 ,9X,' 20    Superpose a curve',/
     5 ,9X,' 88    HELP',/
     2 ,3X,60(1H*),/)

C      IF(.NOT. OKOPN) CALL OPNWRN(1)

      WRITE(6,101) ' * Option  number : '
 101  FORMAT(A20,$)
      READ(5,201,ERR=921) IOPT
 201  FORMAT(I2)
      IF(IOPT.EQ.88) GOTO 88
      GOTO ( 1,921, 3, 4,921, 6,921, 8,99,10,
     > 11,12,13,14,15,16,17,18,19,20 ) IOPT  
      GOTO 921

 1    CONTINUE    
        IZERO=0
        CALL OPNMN(IZERO,
     >               NL,OKOPN,CHANGE,NOMFIC)
        IF(CHANGE) THEN
          IF(KY .EQ. 28) OKBIN=.FALSE.
        ENDIF
      GOTO 921

 3    CONTINUE                      
        CALL PLTOPT(
     >              LM,OKECH,OKBIN,KOPT)
      GOTO 921

 4    CONTINUE
        OKECH=INPECH()
      GOTO 921
        
 6    CONTINUE 
         CALL TRACSY(2,*920)
      GOTO 921

 8    CONTINUE
C        CALL MENVCF
        CALL SAVPLT(*920)
      GOTO 921

 10   CONTINUE 
         CALL PLTEMI(NL,-1)
      GOTO 920

 11   CONTINUE 
         CALL ISOFLD(OKECH,OKOPN,OKVAR,NL,*921)
      GOTO 920

 12   CONTINUE
        CALL CLSCR
        OKECH=.FALSE.
      GOTO 921

 13   CONTINUE 
         CALL SPCTRA(NLOG,NL,LM,OKOPN,CHANGE)
         GOTO 921
C 98      WRITE(6,*) ' ERREUR OPEN ',NOMFIC
C      GOTO 920

 14   CONTINUE 
         CALL LIPS(NLOG,NL,LM,OKOPN,NOMFIC,CHANGE)
      GOTO 921

 15   CONTINUE
      WRITE(6,*) 
      WRITE(6,*) ' KILL/CREATE GRAPHIC SUB-PROCESS (0/1):'
      READ(5,*,ERR=15) KC
      IF(KC.EQ. 0) THEN
        CALL ENDVCF
      ELSEIF(KC.EQ. 1) THEN
        CALL BEGVCF
      ENDIF   
      GOTO 921

 16   CONTINUE 
        CALL SRMEN(NLOG,NL,LM,OKOPN,NOMFIC)
      GOTO 921

 17   CONTINUE 
        CALL CLORBI(NLOG,OKECH,*920)
      GOTO 921

 18   CONTINUE 
        CALL FMATCH(NLOG,NOMFIC)
      GOTO 921

 19   CONTINUE 
        CALL HMATCH(NLOG,NOMFIC)
      GOTO 921

 20   CONTINUE 
C Superpose a curve
        CALL SUPERP(OKECH,NLOG,*920)
      GOTO 921

 88   CONTINUE 
        CALL HELP('1/8')
      GOTO 920

 99   CONTINUE
C      CALL CLMITV       
      IF(OKOPN) THEN
        CLOSE(NL) 
        OKOPN = .FALSE.
      ENDIF
      RETURN 
      END
