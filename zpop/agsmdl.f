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
CDECK GRAPH1 
      SUBROUTINE agsmdl(NLOG,LM,NOMFIC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) NOMFIC
      
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/ OKECH, OKVAR, OKBIN
      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*7, KPOL(2)*9, KDIM(MXVAR)*7
      COMMON/INPVR/ KVAR, KPOL, KDIM
      COMMON/LUN/NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/REBELO/ NRBLT,IPASS,KWRI,NNDES,STDVM
      INCLUDE 'MAXNTR.H'
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      LOGICAL OKOPN, INPECH, CHANGE, KLIPS, KHIST, OK24, OK23
      CHARACTER(14) TXTL
      LOGICAL FIRST

      CHARACTER(1) KLET

      CHARACTER RLIPS,RLIPS0,RLIPS2,RSNPT,RSNPT0,RSNPT2
      CHARACTER RHIST,RHIST0,RHIST2

      SAVE FIRST
      SAVE NL
      SAVE NPASS

      LOGICAL OKREW
     
      LOGICAL IDLUNI
      CHARACTER FNAM*11
      LOGICAL FIRST2
      SAVE FNAM, NLIPS, FIRST2

      PARAMETER (ITWO=2)
     
      integer debstr, finstr
      character(400) cmmnd

      DATA  OKOPN / .FALSE. /
      DATA CHANGE / .TRUE./
      DATA FIRST / .TRUE. /
      DATA RSNPT,RSNPT2 / 'N','Y' /
      DATA RLIPS,RLIPS2 / 'Y','N' /
      DATA RHIST,RHIST2 / 'Y','N' /
      DATA OKREW / .TRUE. /

      DATA FNAM /'zpop.lips'/
      DATA FIRST2 /.TRUE./

c       nlips = 77
c         goto 921
c         write(*,*) 'Pgm graph7 ' 
c          read(*,*)

c      IF(FIRST2) THEN
c        IF (IDLUNI(NLIPS)) THEN
c          OPEN(UNIT=NLIPS,FILE=FNAM,ERR=997)
c          CLOSE(UNIT=NLIPS,STATUS='DELETE')
c          OPEN(UNIT=NLIPS,FILE=FNAM,STATUS='NEW',ERR=997)
c          WRITE(6,*)
c     >      ' Pgm graph7. Opened ',FNAM(debstr(FNAM):finstr(FNAM))
c        ELSE
c          WRITE(6,*) ' *** Problem in graph7 : No idle unit number ! '
c          WRITE(6,*) '     No idle unit number !  Forced to 77 '
c          NLIPS=77
c        ENDIF
c        FIRST2 = .FALSE.
c      ENDIF

      GOTO 921

 920  CONTINUE      
      CALL FBGTXT
      I=IDLG('('' Press RETURN for more :'')','    ',1)

 921  CONTINUE
      CALL FBGTXT
      CALL HOMCLR

      WRITE(6,100)     !!!!NOMFIC,KVAR(KX),KX,KVAR(KY),KY
 100  FORMAT(//,3X,60('*'),//,20X,' MENU - ZgoubiFromSnaprampCmd :' ,//
     > ,9X,'        ' ,/
c     1 ,9X,'  1    OPEN  FILE - current is ',A ,/
c     2 ,9X,'  2    VARIABLES  TO  PLOT ' ,/
c     2 ,9X,'            PRESENT:  H-AXIS :',A,'(KX=',I2,')',/
c     2 ,9X,'                      V-AXIS :',A,'(KY=',I2,')',/  
c     3 ,9X,'  3    PLOT  OPTIONS',/
c     4 ,9X,'  4    MANUAL  SCALES ',/
c     5 ,9X,'  5    AUTOMATIC  SCALES ',/
     7 ,9X,'  7    ZgoubiFromSnaprampCmd   ',/
c     8 ,9X,'  8    Save  graphic  screen ',/
     9 ,9X,'  9    EXIT  THIS  MENU     ',/
c     > ,9X,' 12    ERASE  DISPLAY    ',/
     > ,9X,' 15    KILL/CREATE  GRAPHIC  SUB-PROCESS      ',/
c     > ,9X,' 20    Superpose a curve',/
c     > ,9X,' 77    **  PLOT  Y V.S. X, turn by turn  **    ',/
     > ,9X,' 88    HELP',/
     > ,3X,60('*'),/)

c      IF(FIRST) THEN
c        CALL PLINIT('////MENU07////',
c     >                               NL,OKOPN,CHANGE,NOMFIC)
c        FIRST = .FALSE.
c      ENDIF
c      IF(.NOT. OKOPN) CALL OPNWRN(1)

      WRITE(6,101) ' * Option  number : '
 101  FORMAT(A20)
      READ(5,201,ERR=921) IOPT
 201  FORMAT(I2)
      IF(IOPT.EQ.88) GOTO 88
      GOTO (921,921,921,921,921,921, 7,921, 99) IOPT  
      GOTO 921

 7    CONTINUE 
      cmmnd = 
     >'/rap/lattice_tools/zgoubi/AgsZgoubiModel/ZgoubiFromSnaprampCmd'
      write(*,fmt='(a)') cmmnd(debstr(cmmnd):finstr(cmmnd))
      call system(cmmnd)
      GOTO 920
    
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

 20   CONTINUE 
C Superpose a curve
        CALL SUPERP(OKECH,NLOG,*920)
      GOTO 920

 88   CONTINUE 
        CALL HELP('1/7')
      GOTO 920

 99   CONTINUE
C      CALL CLMITV       
      IF(OKOPN) THEN
        CLOSE(NL) 
        OKOPN = .FALSE.
      ENDIF
      RETURN 

 997  WRITE(6,*) ' Error upon OPEN ',FNAM
      RETURN 

      END
