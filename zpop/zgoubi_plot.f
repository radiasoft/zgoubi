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

      CHARACTER * 9   DMY,HMSF  
      CHARACTER * 80  NOMFIC
      LOGICAL FIRST
      SAVE FIRST
                                 
      DATA FIRST / .TRUE. /

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

 998     I=I+1
         OPEN(UNIT=NLOG,FILE= 'zpop.log' ,STATUS='OLD')
         CLOSE(UNIT=NLOG,STATUS='DELETE')
         IF(I.LT. 10) GOTO 898
         WRITE(6,*) ' SBR ZGOUBI : ERREUR  OPEN FICHIER  zpop.log'

         CONTINUE
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

 7    CONTINUE                      
         CALL GRAPH7(NLOG, LM, NOMFIC)
         GOTO 21

 8    CONTINUE                      
         CALL GRAPH8(NLOG, LM, NOMFIC)
         GOTO 21

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
