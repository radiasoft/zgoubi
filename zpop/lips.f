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
CDECK LIPS
      SUBROUTINE LIPS(NLOG,NL,LM,OKOPN,NOMFIC,CHANGE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKOPN,CHANGE
      CHARACTER*(*) NOMFIC
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN
      COMMON/LUN/NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      INCLUDE 'MAXNTR.H'
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR

      DIMENSION YM(3), YPM(3), U(3), A(3), B(3)
      DIMENSION YMX(6), YPMX(6)
      CHARACTER REP
      DATA KPS / 1 /
      LOGICAL BINARY, BINARF
      INCLUDE 'FILFAI.H'

      DIMENSION XSIGU(3)
      DATA XSIGU / 3*1.D0 /
      SAVE XSIGU

      MODSTO = 2

      IF(.NOT.OKOPN)
     > CALL OPNDEF(NFAI,FILFAI,NL,
     >                            NOMFIC,OKOPN) 

      NPTR=NPTS
      GOTO 21

 20   CONTINUE      
      CALL FBGTXT
      WRITE(*,FMT='(/,''  Press RETURN for more'')') 
      READ(*,200,ERR=20) REP 
 200  FORMAT(A1)

 21   CONTINUE
      CALL FBGTXT
      CALL HOMCLR

      CALL READC5(NT1,NT2)
      WRITE(*,104) NOMFIC,NT1,NT2,LM,NPTR,NPTS
 104  FORMAT(/,5X,' MENU  -  fit of phase-space ellipses : ',//
     1,/,5X,' 1  OPEN  FILE     -  current is ',A
     2,/,5X,' 2  Particle  #    (',I6,' to ',I6,')'
     >,/,5X,'    and, # of the recording element (-1 for all)  (',I6,')'
     3,/,5X,' 3  Number  of  phase-space  points  to  read   (',I6,')'
     >,/,5X,'                                and  to  fit    (',I6,')'
     >)
C      IF(NT.EQ. -1) WRITE(*,105) KPS
      IF(NT2.GT.NT1) WRITE(*,105) KPS
 105  FORMAT(5X,' 5  Fit initial or final phase-space coordinates',
     5' (0/1)  (now ',I1,') : ',$)
      WRITE(*,106) 
 106  FORMAT(/,
     6   5X,' 6  **  Compute fitting ellipse parameters   ** '                           
     7,/,5X,' 7  **  Plot ellipses  ** '                           
     8,/,5X,' 8  Print screen'
     9,/,5X,' 9  EXIT  THIS  MENU '
     X,/,5X,'12  ERASE  DISPLAY    '
     >,/)

      WRITE(*,100)
 100  FORMAT('$  Option  number : ',$) 
      READ(*,108,ERR=21) IOPT
108   FORMAT(I2)

      BINARY=BINARF(NL)

      IF(.NOT. OKOPN .AND. (IOPT.EQ.6 .OR. IOPT.EQ.7)) THEN
        CALL OPNWRN(1)           
        GOTO 20
      ELSE
        GOTO ( 1, 2, 3,21,5, 6, 7, 8, 9,21,21,12) IOPT  
        GOTO 21
      ENDIF

 1    CONTINUE                      
        CALL OPNMNL(NFAI,
     >                  NL,NOMFIC,OKOPN,CHANGE)
      GOTO 20

 2    CONTINUE
      NT10 = NT1
      NT20 = NT2
      CALL READC6(5)
      CALL READC4(5)
      CALL READC3(KL1,KL2)
      LM = KL1

      IF((NT1 .NE. NT10  .OR. NT2 .NE. NT20) .OR. LM .NE. LM0) 
     >                                           CHANGE = .TRUE.
      GOTO 21

 3    CONTINUE    
      NPTR0=NPTR
      WRITE(*,130) NTRMAX, NPTR
 130  FORMAT
     >(/,' Read phase-space coordinates up to (max=',I7,')   (now ',
     > I7,') : ',$)
      READ(5,*,ERR=3) NPTR
      IF(NPTR.LE. 0 .OR. NPTR.GT. NTRMAX) GOTO 3
      IF(NPTR.GT. NPTR0) CHANGE=.TRUE.
      IF(NPTS.GT. NPTR) NPTS=NPTR

 31   WRITE(*,140) NPTS
 140  FORMAT(/,'  Number of phase-space coordinates'
     >,1X,'to fit (max =  # read)   (now ',I6,') : ',$)
      READ(5,*,ERR=31) NPTS
      IF(NPTS.LE. 0) GOTO 31
      IF(NPTS.GT. NPTR) NPTS=NPTR
      GOTO 21

 5    CONTINUE
      KPS0=KPS
      WRITE(*,150) KPS
 150  FORMAT(/,'  Initial or final coordinates (0/1) (now ',I1,') : ',$)
      READ(5,*,ERR=5) KPS
      IF(KPS.NE.0 .AND. KPS.NE.1) GOTO 5 
      IF(KPS.NE.KPS0) CHANGE=.TRUE.
      GOTO 21

 6    CONTINUE
          IF(CHANGE) THEN
            CALL STORCO(MODSTO,NL,KPS,
     >                                          NPASS)
            CHANGE=.FALSE.
          ENDIF
          IF(NPTR .GT. 0) THEN
            IF(NPTS.GT. NPTR) NPTS=NPTR
            CALL LPSFIT(NLOG,1,LM,
     >                            YM,YPM,YMX,YPMX,U,A,B,*60,*60)
 60         CONTINUE
          ENDIF
          CALL LPSCNT(YM,YPM,U,A,B,XSIGU,NLOG,
     >                                        NCOUNT)
      GOTO 20

 7    CONTINUE
        CALL LPSDRW(NLOG,YM,YPM,YMX,YPMX,U,A,B,KPS,XSIGU,*21)
      GOTO 20

 8    CONTINUE
C        CALL MENVCF
        CALL SAVPLT
      GOTO 21

 9    RETURN

 12   CONTINUE
        CALL CLSCR
        OKECH=.FALSE.
      GOTO 21

C 98   WRITE(*,*) ' ERREUR OPEN ',NOMFIC
C      RETURN  
      END
