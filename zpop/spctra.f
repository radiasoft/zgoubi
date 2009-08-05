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
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE SPCTRA(NLOG,NL,LM,OKOPN,CHANGE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKOPN, CHANGE
C----- PLOT SPECTRUM     
      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      PARAMETER (NCANAL=2500)
      COMMON/SPEDF/BORNE(6),SPEC(NCANAL,3),PMAX(3),NC0(3)
      INCLUDE 'MAXNTR.H'
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR

      DIMENSION YM(3), YPM(3), U(3), A(3), B(3), YNU(3)
      DIMENSION YMX(6), YPMX(6)
 
      LOGICAL OKECH
      CHARACTER REP, NOMFIC*80
      LOGICAL BINARY, BINARF
      CHARACTER HVL(3)*12

      SAVE NT
      LOGICAL OPN

      LOGICAL IDLUNI, OKKT5

      INCLUDE 'FILFAI.H'

      SAVE NPASS

      DATA NT / -1 /

      DATA OPN / .FALSE. /
      DATA HVL / 'Horizontal', 'Vertical', 'Longitudinal' /
    
      IF(.NOT.OKOPN) 
     > CALL OPNDEF(NFAI,FILFAI,NL,
     >                            NOMFIC,OKOPN) 

      MODSTO = 1
      NPTR=NPTS

      OKECH = .FALSE.

      GOTO 21

 20   CONTINUE      
      CALL FBGTXT
      WRITE(*,FMT='(/,''  Press RETURN for more'')') 
      READ(*,200,ERR=20) REP 
 200  FORMAT(A1)

 21   CONTINUE
      CALL FBGTXT
      CALL HOMCLR

      IF(NT.EQ.-1) THEN
        CALL READC6B(1,NPTS)
      ELSE
        CALL READC6B(NT,NT)        
      ENDIF

      WRITE(*,104) NOMFIC,MXT,NT,LM,NPTR,NPTS
 104  FORMAT(5X,' MENU  -  Fourier  spectrum : ',/
     1,/,5X,' 1  OPEN  FILE     -  current is ',A
     2,/,5X,' 2  Particle  # (1 - ',I6,',   -1 for all)       (',I6,')'
     >,/,5X,'    and, # of the recording element (-1 for all)  (',I6,')'
     3,/,5X,' 3  Number  of  phase-space  points  to  read  (',I6,')'
     >,/,5X,'          and number  of  points  to  analyse  (',I6,')'
     5,/,5X,' 5  Define Nu-min / Nu-max   '                
     5,/,5X,'          and  number  of  channels   '
     6,/,5X,' 6  **  Compute  tunes, beta, etc.  ** '                          
     7,/,5X,' 7  **  Plot  ** '                          
     8,/,5X,' 8  Print screen'
     9,/,5X,' 9  EXIT  THIS  MENU '
     2,/,5X,'12  ERASE  DISPLAY    '
     >,/)

      WRITE(*,100)
 100  FORMAT('$  Option  number : ')
      READ(*,108,ERR=21) IOPT
108   FORMAT(I2)
      BINARY=BINARF(NL)
      IF(.NOT. OKOPN .AND. (IOPT.EQ.6 .OR. IOPT.EQ.7)) THEN
        CALL OPNWRN(1)           
        GOTO 20
      ELSE
        GOTO ( 1, 2, 3,21, 5, 6, 7, 8, 9,21,21,12) IOPT  
        GOTO 21
      ENDIF

 1    CONTINUE    
        CALL OPNMNL(NFAI,
     >                   NL,NOMFIC,OKOPN,CHANGE)
      GOTO 20

 2    CONTINUE
      NT0 = NT
 23   WRITE(*,FMT=
     > '(''  # of the particle to be analized ( 1-'',I6,'' ): '')') MXT 
      WRITE(*,FMT=
     > '(''      ( Default = increment )'',I3)') NT
      READ(5,FMT='(I4)',ERR=23) NT
      IF(NT .GT.  0) THEN 
        IF(NT.LE.MXT) THEN
          CALL READC6B(NT,NT)
        ELSE
          GOTO 23
        ENDIF
      ELSEIF(NT .EQ.  0) THEN 
        IF(NT0.EQ.-1) THEN
          NT0 = NT
          GOTO 23
        ELSE
          NT = NT0 + 1
          CALL READC6B(NT,NT)
        ENDIF
      ELSEIF(NT.LT. 0) THEN 
        IF(NT.EQ.-1) THEN 
          CALL READC6B(1,NPTS)
        ELSEIF(NT.LT.-1) THEN
          GOTO 23
        ENDIF
      ENDIF
      CALL READC4(5)
      CALL READC3(KL1,KL2)
      LM = KL1

      IF(NT .NE. NT0 .OR. LM .NE. LM0) CHANGE = .TRUE.
      GOTO 21

 3    CONTINUE    
      NPTR0=NPTR
      WRITE(*,130) '  NUMBER OF PHASE-SPACE COORDINATES TO READ',
     >  '  (max is ',NTRMAX,' ;  now ',NPTR,')  :'
 130  FORMAT(/,A43,A9,I6,A8,I6,A4) 

      READ(5,*,ERR=3) NPTR
      IF(NPTR.LE. 0 .OR. NPTR.GT. NTRMAX) GOTO 3
      IF(NPTR.GT. NPTR0) CHANGE=.TRUE.
      IF(NPTS.GT. NPTR) NPTS=NPTR

 31   WRITE(*,131) NPTS
 131  FORMAT(/,'  NUMBER OF PHASE-SPACE COORDINATES'
     >,1X,'TO ANALYSE (max =  # read):',I7)
      READ(5,*,ERR=31) NPTS
      IF(NPTS.LE. 0) NPTS=NPTR
      IF(NPTS.GT. NPTR) NPTS=NPTR
      GOTO 21

 5    CONTINUE

      DO 52 INU = 1, 5, 2

        JNU = 1 + INU/2
        ANUI = BORNE(INU)
        ANUF = BORNE(INU+1)
        WRITE(*,FMT='(A,''  motion :'')') HVL(JNU)
 53     WRITE(*,150) ANUI, ANUF
 150    FORMAT(/,' Frequency spectrum range ;  Nu-min, max :',2F8.5)
        READ(5,*,ERR=53) ANUI, ANUF
        IF (ANUF.LE. ANUI) GO TO 53
        BORNE(INU)=ANUI
        BORNE(INU+1)=ANUF

 51     WRITE(*,151) NC0(JNU)
 151    FORMAT(1X,' NB OF CHANNELS (<2500): ',I5)
        READ(5,*,ERR=51) NC0(JNU)
        IF (NC0(JNU) .LE. 0 .OR. NC0(JNU) .GT. NCANAL) GO TO 51
        OKECH = .FALSE.
 52   CONTINUE
      GOTO 21

 6    CONTINUE
        IF(.NOT. OPN) THEN
          IF(NT.EQ.-1) THEN
            IF (IDLUNI(IUN)) THEN
              OPEN(UNIT=IUN,FILE='zpop.outLips',ERR=699)
              OPN = .TRUE.
            ELSE
              GOTO 698
            ENDIF
          ENDIF
        ENDIF

          IF(NT.EQ.-1) THEN
            KPR = 2
            KT = 1
            CALL READC6B(KT,KT)
          ELSE
            KPR = 1
          ENDIF
         
 62       CONTINUE
          IF(CHANGE) THEN
            NPTR = NTRMAX
            NPTS=NPTR
            CALL STORCO(MODSTO,NL,1  ,
     >                                          NPASS)
c            write(*,*) ' spctra  npass : ', npass

            CHANGE=.FALSE.
            IF(NPTR.GT.0) THEN
              IF(NPTS.GT. NPTR) NPTS=NPTR
            ELSE
              NPTR = NPTS
              IF(NT.EQ.-1) KT = 1 
             GOTO 69
            ENDIF
          ENDIF
          IF(NPTR .GT. 0) THEN
            CALL LPSFIT(NLOG,KPR,LM,
     >                              YM,YPM,YMX,YPMX,U,A,B,*60,*60)
 60         CONTINUE
            WRITE(*,FMT='(/,A,/)') '  Busy, computing tunes...'
            CALL SPEANA(YM,BORNE,NC0,
     >                               YNU,SPEC,PMAX)
            CALL SPEPR(NLOG,KPR,LM,NT,NPTS,YM,YPM,YNU,PMAX,NC0)
            IF(OKKT5(KT)) THEN
              IF(KPR.EQ.2) CALL SPEIMP(IUN,YNU,BORNE,U,KT)
            ENDIF
          ENDIF

          IF(NT.EQ.-1) THEN
            KT = KT+1
            CALL READC6B(KT,KT)
            CHANGE = .TRUE.
            GOTO 62
          ENDIF

 69       CONTINUE
          IF(OPN) THEN
            WRITE(IUN,*) ' XXNU,ZZNU, (U(I),I=1,3),COOR(KT,6), KT, NPTS'
            CALL FLUSH2(IUN,.FALSE.)
            CHANGE = .TRUE.
          ENDIF
      GOTO 20
 698  WRITE(6,*) ' *** Problem : No idle unit for zpop.out_Lips '
      GOTO 20
 699  WRITE(6,*) ' *** Problem at OPEN zpop.out_Lips '
      GOTO 20

 7    CONTINUE
          CALL SPEDRW(NT,BORNE,YNU,PMAX,SPEC,NC0,OKECH)
          GOTO 21

 8    CONTINUE
C        CALL MENVCF
        CALL SAVPLT
      GOTO 21

 9    RETURN

 12   CONTINUE
        CALL CLSCR
      GOTO 21

      END
