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
C  -------
      PROGRAM ZGOUBI_PLOT
C PGM PRINCIPAL DU GRAPHIQUE PACKAGE DESTINE AU TRAITEMENT
C DES OUTPUT DE ZGOUBI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN    
      COMMON/LUN/                        L1  ,L2  ,L3  ,L4  ,L5  ,L6  

      CHARACTER(10)   DMY
      CHARACTER(9)   HMSF  
      CHARACTER(80)  NOMFIC
      LOGICAL FIRST, FIRST1, EXS, OK, IDLUNI
      SAVE FIRST, FIRST1 
      CHARACTER(400) WRKDIR
      INTEGER DEBSTR, FINSTR
                                      
      DATA FIRST, FIRST1 / .TRUE. , .TRUE. /

      CALL SYSTEM('rm -f echoDir ; pwd | cat > echoDir ')
      OK = IDLUNI(
     >            LDIR)
      OPEN(LDIR,FILE='echoDir')
      READ(LDIR,FMT='(A)') WRKDIR
      CLOSE(LDIR)

      WRITE(*,*) 'Pgm zgoubi_plot.  Working directory : ', 
     >   WRKDIR(DEBSTR(WRKDIR):FINSTR(WRKDIR))

      CALL INIGR(
     >           LM, NOMFIC)
      CALL INIGR1(
     >            NLOG)
      
      L1 = NDAT
      L2 = NRES
      L3 = NPLT
      L4 = NFAI
      L5 = NMAP
      L6 = NSPN

      I=0
      INQUIRE(FILE='zpop.log',EXIST=EXS)
 898  continue
      IF(EXS) THEN
        OPEN(UNIT=NLOG,FILE='zpop.log',STATUS='OLD',ERR=998
     >  ,IOSTAT=IOS)
        CLOSE(NLOG,STATUS='DELETE')
      ENDIF
c         write(*,*) 'Pgm zgoubi_plot : exs=' ,exs
c         read(*,*)
      OPEN(UNIT=NLOG,FILE='zpop.log',ERR=998,IOSTAT=IOS)
      
      GOTO 21
 
 20   CONTINUE      
      I=IDLG('('' Press RETURN for more :'')','    ',1)

 21   CONTINUE
      IF(FIRST) THEN
        CALL BEGVCF
        FIRST = .FALSE.
      ENDIF
      CALL FBGTXT
      CALL MNZGRA(IOPT,WRKDIR) 
      GOTO (21,21,21,21,21,21, 7, 8,21,99,11,12) IOPT  
      GOTO 21

 998     continue

         CLOSE(NDAT) 
         CLOSE(NRES) 
         CLOSE(NFAI) 
         CLOSE(NPLT) 
         CLOSE(NSPN)
         CLOSE(NMAP) 

         CALL DATE2(DMY)                     
         CALL TIME2(HMSF)                   
         
         IF(.NOT. FIRST1) THEN
           WRITE(6,100) DMY,HMSF
 100       FORMAT(/,'  Job  ended  on  ',A,',  at  ',A,/)
           GOTO 20      
         else
           first1 = .false.
           goto 21
         endif

 7    CONTINUE                      
c         write(*,*) 'Pgm zgoubi_plot : call graph7' 
c         read(*,*)
         CALL GRAPH7(NLOG, LM, NOMFIC,wrkDir)
         GOTO 21

 8    CONTINUE                      
         CALL GRAPH8(NLOG, LM, NOMFIC,WRKDIR)
         GOTO 21

 11      CONTINUE                      
         CALL AGSMDL(NLOG, LM, NOMFIC)
         GOTO 21

 12      CONTINUE                      
         CALL AGSDDQ(NLOG, LM, NOMFIC)
         GOTO 21

99    CONTINUE
         CALL ENDVCF
 

         STOP  
      END 
