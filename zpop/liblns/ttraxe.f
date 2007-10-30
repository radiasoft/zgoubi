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
C*****TEST TRAXE
      CHARACTER*81 TAMP
      CHARACTER*20 PRN
      INTEGER*2 ISYM1(10),ISYM2(12)
      COMMON/EDIT/TAMP
      DATA ISYM1/4,-1000,-577,1000,-577,0,1155,-1000,-577,0/
      DATA ISYM2/5,-1,-1,1,-1,1,1,-1,1,-1,-1,0/
      DATA MODE,X1,X2,Y1,Y2,NPOINT,NTYPE/2,-4.,4.,-1.,1.,100,1/
      CALL TXTFBG
      CALL INIVCF
C     CALL INISCA(1535,1139)
C      CALL INISCA(1024, 812)
C      CALL INISCA(880,670)
C      CALL INISCA(512, 384)
      CALL RPSUPP(.TRUE.)
      I=IFMES(0)
1     CONTINUE
      CALL FBGTXT
C      CALL SCRCLR(1,1)      
      WRITE(6,100)
  100 FORMAT(' ENTREZ MODE,X1,X2,Y1,Y2,NP,TYP OU -1 POUR TERMINER:')
      WRITE(TAMP,'(1X,I2,4G13.4,I5,I3)')MODE,X1,X2,Y1,Y2,NPOINT,NTYPE
      CALL RPSLDS(1)
      CALL RPSLDI(MODE)
      IF(MODE.EQ.-1) THEN
	 CALL FINVCF
	 STOP 1
      ENDIF
      CALL RPSLDE(X1)
      CALL RPSLDE(X2)
      CALL RPSLDE(Y1)
      CALL RPSLDE(Y2)
      CALL RPSLDI(NPOINT)
      CALL RPSLDI(NTYPE)
      CALL TXTFBG
C      CALL TRSTRT
      CALL VECTF(0,0,0)
C      CALL TRSTRT
      CALL TRAXE(X1,X2,Y1,Y2,MODE)
C      CALL TRSTRT
      CALL INIVPH(NTYPE)
      IF(NTYPE.GT.10) THEN
         CALL DEFMKR(10,NTYPE-10)
         CALL INIVPH(9)
      ENDIF
      DELX=(X2-X1)/NPOINT
      Z=SIN(X1)
      CALL VECTPH(X1,Z,4)
      NPNT=NPOINT+1
      DO 10 I=1,NPNT
        IF(IFMES(3).EQ.1) THEN
          CALL FBGTXT
          WRITE(*,*) "Interruption par l'opérateur"
          GOTO 11
        ENDIF 
      DO II=1,10000
	XX=SIN(FLOAT(I))
      END DO
      T=X1+(I-1)*DELX
      Z=SIN(T)
      CALL VECTPH(T,Z,2)
      IF(NTYPE.NE.9)GOTO 10
      CALL INIVPH(1)
      COEF1=0.020
      COEF2=0.06*(1.+ABS(Z)/0.5)
      CALL TRSYMB(T,Z,2,ISYM1,0,0,COEF1,COEF1)
      CALL TRSYMB(T,Z,2,ISYM2,0,1,16.,COEF2)
      CALL INIVPH(NTYPE)
   10 CONTINUE
   11 CONTINUE
      X=X1+0.1*(X2-X1)
      Y=Y1+0.9*(Y2-Y1)
      CALL DEFCAR(4,0,0)
      CALL TRTEXT(X,Y,'Y=SIN(X)',1)
      CALL TRTEXT(0.,0.,
     *'12345678901234567890123456789012345678901234567890123456789012345
     *678901234',0)
      CALL DEFCAR(3,0,0)
      CALL TRTEXT(0.,20.,
     *'12345678901234567890123456789012345678901234567890123456789012345
     *678901234567890',0)
      CALL DEFCAR(2,0,0)
      CALL TRTEXT(0.,40.,
     *'12345678901234567890123456789012345678901234567890123456789012345
     *678901234567890',0)
      CALL DEFCAR(1,0,0)
      CALL TRTEXT(0.,60.,
     *'12345678901234567890123456789012345678901234567890123456789012345
     *678901234567890',0)
C      CALL TRSTOP

      CALL DEFCAR(2,1,0)
      CALL FBGTXT
      CALL TRGETD(X1,X2,DIVX,NTX)
      CALL TRGETD(Y1,Y2,DIVY,NTY)
      WRITE(6,*) ' DIVX=',DIVX,' DIVY=',DIVY,' NTX=',NTX,' NTY=',NTY
      CALL DLGDEF(.TRUE.)
      I=IDLG('('' COPIE ECRAN (F=fichier,P=impr.,[N]=non ) :'')'
     *         ,'N   F   P   ',3)
      IF (I.EQ.2) THEN
	  CALL SAVECR('test.ps')
      ELSE IF(I.EQ.3) THEN
          TAMP=' NOM DE L''IMPRIMANTE:'
          CALL RPSUPP(.FALSE.)
          CALL RPSLDS(0)
          CALL RPSLDT(PRN)
          CALL RPSUPP(.TRUE.)
          CALL SAVECR(PRN//':')
      ENDIF
      GOTO 1
      END
