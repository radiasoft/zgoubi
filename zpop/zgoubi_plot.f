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
      PROGRAM ZGOUBI_PLOT
C PGM PRINCIPAL DU GRAPHIQUE PACKAGE DESTINE AU TRAITEMENT
C DES OUTPUT DE ZGOUBI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN    
      COMMON/LUN/                        L1  ,L2  ,L3  ,L4  ,L5  ,L6  

      CHARACTER * 9   DMY,HMSI,HMSF  
      CHARACTER * 80  CMND, NOMFIC
      LOGICAL FIRST
      SAVE FIRST
                                 
      CHARACTER * 27 COMMND

      DATA FIRST / .TRUE. /

      CALL BLOCK
      CALL INIGR(
     >           NLOG, LM, NOMFIC)
      
      L1 = NDAT
      L2 = NRES
      L3 = NPLT
      L4 = NFAI
      L5 = NMAP
      L6 = NSPN

      I=0
 898  OPEN(UNIT=NLOG,FILE='zpop.log',STATUS='NEW',ERR=998
     >,IOSTAT=IOS)
      IF(IOS.NE.0) THEN 
        CLOSE(NLOG,STATUS='DELETE')
        GOTO 898
      ENDIF

      GOTO 21
 
 20   CONTINUE      
      I=IDLG('('' Press RETURN for more :'')','    ',1)

 21   CONTINUE
      IF(FIRST) THEN
        CALL BEGVCF
        FIRST = .FALSE.
      ENDIF
      CALL FBGTXT
      CALL MNZGRA(IOPT) 
      GOTO (21,21,21,21,21,21, 7, 8,21,99,21) IOPT  
      GOTO 21

 1    CONTINUE

         OPEN(UNIT=NDAT,FILE='zgoubi.dat' ,STATUS='OLD',ERR=996)
         I=0
 897     OPEN(UNIT=NRES,FILE='zgoubi.res' ,STATUS='NEW',ERR=997
     >   ,IOSTAT=IOS)
         IF(IOS.NE.0) THEN 
           CLOSE(NRES,STATUS='DELETE')
           GOTO 897
         ENDIF

         CALL DATE2(DMY)                     
         CALL TIME2(HMSI) 
         WRITE(NRES,103) DMY,HMSI 
         WRITE(6   ,103) DMY,HMSI
 103     FORMAT(/,'  Job  started  on  ',A9,',  at  ',A9,/)

C         CALL ZGOUBI(200,.TRUE.)
         COMMND = '~/meot/zgoubi/source/zgoubi'
         CALL SYSTEM(COMMND)
         GOTO 10

 996     WRITE(6,*) ' SBR ZGOUBI : ERREUR  OPEN FICHIER  zgoubi.dat '
         GOTO 10
 997     I=I+1
         OPEN(UNIT=NRES,FILE= 'zgoubi.res' ,STATUS='OLD')
         CLOSE(UNIT=NRES,STATUS='DELETE')
         IF(I.LT. 10) GOTO 897
         WRITE(6,*) ' SBR ZGOUBI : ERREUR  OPEN FICHIER  zgoubi.res '
         GOTO 10
 998     I=I+1
         OPEN(UNIT=NLOG,FILE= 'zpop.log' ,STATUS='OLD')
         CLOSE(UNIT=NLOG,STATUS='DELETE')
         IF(I.LT. 10) GOTO 898
         WRITE(6,*) ' SBR ZGOUBI : ERREUR  OPEN FICHIER  zpop.log'

 10      CONTINUE
         CLOSE(NDAT) 
         CLOSE(NRES) 
         CLOSE(NFAI) 
         CLOSE(NPLT) 
         CLOSE(NSPN)
         CLOSE(NMAP) 

         CALL DATE2(DMY)                     
         CALL TIME2(HMSF)                   
         WRITE(6,100) DMY,HMSF
 100     FORMAT(/,'  Job  ended  on  ',A9,',  at  ',A9,/)
         GOTO 20      

 2    CONTINUE            
         CMND='BATCH ZGOUBI'
         IOUPI=LIBSPA(CMND) 
         GOTO 20 

 3    CONTINUE
         CALL MNFILE          
         GOTO 21

 7    CONTINUE                      
         CALL GRAPH1(NLOG, LM, NOMFIC)
         GOTO 21

 8    CONTINUE                      
         CALL GRAPH2(NLOG, LM, NOMFIC)
         GOTO 21

 11   CONTINUE
         CMND= '   SHOW  QUEUE  ' 
         IOUPI=LIBSPA(CMND)
      GOTO 20 

99    CONTINUE
         CALL ENDVCF
C         WRITE(6,102) 
C 102     FORMAT('$ SCRATCH  zgoubi.plt .  CONFIRM (N/O) : ')
C         READ(5,FMT='(A1)',ERR=99) REP
C         IF    (REP.EQ. 'O' .OR. REP.EQ. 'o') THEN
C            CMND= 'DELETE/NOCONFIRM zgoubi.plt;*'
C            IOUPI=LIBSPA(CMND)
C         ELSEIF(REP.EQ. 'N' .OR. REP.EQ. 'n') THEN
C            WRITE(6,*) '    Well ...'
C         ENDIF 
 
         STOP  
      END 
