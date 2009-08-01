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
      SUBROUTINE PLTOPT(
     >                  LM,OKECH, OKBIN,IO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKECH, OKBIN

      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*7, KPOL(2)*9, KDIM(MXVAR)*7
      COMMON/INPVR/ KVAR, KPOL, KDIM
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      PARAMETER (MXB=1000)
      CHARACTER KAX(2)*14, KLIS(2)*3

      LOGICAL  OKXAV, OKYAV
      SAVE   OKXAV, OKYAV
      LOGICAL  OKX12, OKY12
      SAVE   OKX12, OKY12

      CHARACTER*1 KLET

      SAVE AX,PX,BX,AY,PY,BY,IBXY

      LOGICAL IDLUNI

      DATA KAX/ 'Full  scale   ',  'Proportionnal ' /
      DATA KLIS / 'No ', 'Yes'/
      DATA OKXAV,OKYAV / .FALSE., .FALSE./
      DATA OKX12,OKY12 / .FALSE., .FALSE./
      DATA AX,PX,BX,AY,PY,BY / 1.D0, 1.D0, 0.D0, 1.D0, 1.D0, 0.D0 /
      DATA IBXY,IUN / 0, 0 /
      
 21   CONTINUE
      CALL HOMCLR

      
      CALL READC1(KP1,KP2)
      CALL READC5(KT1,KT2)
      CALL READC9(KKEX,KLET)
      WRITE(6,116) LIS,IAX,NB,KKEX,KLET,KP1,KP2,LM,KT1,KT2,OKXAV,OKYAV,
     > OKX12,OKY12
 116  FORMAT(/,5X,' PLOT OPTIONS: '
     >,//,'   PLOT OPTIONS:'
     >, /,8X,' 2: List X & Y ( No <-> Yes     - FLIP-FLOP (',I1,') )'
     >, /,8X,' 3: Axes  ( Full <-> Proportional - FLIP-FLOP (',I1,') )'
     >, /,8X,' 4: Histogram with ',I4,'  bins '
     >, /,8X,' 5: Select particle category * KEX/KTag =',I3,'/',A1
     >, /,8X,' 6: Select  pass  #.  Now : ',I6,'  -  ',I6
     >, /,8X,' 7: Select  optical  element  #  (-1 = all).  Now : ',I5
     >, /,8X,' 8: Select  trajectory  #.  Now : ',I6,'  -  ',I6
     >, /,8X,' 9: EXIT' 
     >, /,8X,'10: Line/marker  style         ' 
     >, /,8X,'11: Reset S coordinate at each pass '
     >, /,8X,'12: Tag positioning coefficient '
     >, /,8X,'13: Change H axis to x-<x> :', L1,' V axis to y-<y> :', L1
     >, /,8X,'14: Change H/V variables :  X -> f(X)  and/or  Y -> f(Y) '
     >, /,8X,'15: Plot dx :', L1,'   dy :', L1
C     >, /,8X,'16: Plot <y> +/- sig_y vs. x (e.g., Env. vs. Ipass) :', L1
     >, /,8X,'26: Path extrapolation in field maps'
     >, /,8X,'88: more...' 
     >,/)

      WRITE(6,100) ' * Option  number : '
 100  FORMAT(A20)
      READ(5,*,ERR=21) IO
      GOTO(21, 2,3 ,4 ,5 ,6 ,7 ,8 ,98,10,11,12,13,14,15) IO
      IF(IO.EQ.26) GOTO 26
      IF(IO.EQ.88) GOTO 88
      GOTO 21

 2    CONTINUE
C------ List/NoList X & Y (2/1)
       IF(LIS .EQ. 1) THEN
          LIS = 2   ! list
        ELSE
          LIS = 1   ! no list
        ENDIF 
      GOTO 98

 3    CONTINUE
C-------- Axes  ( Full <-> Proportional)
        OKECH = .FALSE.
        IF(IAX .EQ. 1) THEN
          IAX=2
          CALL TRAXPI
        ELSE
          IAX = 1
          CALL TRAXPJ
        ENDIF 
      GOTO 98

 4    CONTINUE
        NB0=NB
 41     WRITE(6,FMT='('' # bins of histo ( MAX='',I4,'') : '')') MXB
        READ(5,*,ERR=41) NB
        IF(NB .GT. MXB) GOTO 41
        IF(NB .NE. NB0) OKBIN=.FALSE.
      GOTO 98

 5    CONTINUE
 51     WRITE(6,FMT='('' KEX value - 99 for any (now  : '',I4,
     >                                          '') : '')') KKEX
        READ(5,*,ERR=51) KKEX
 52     WRITE(6,FMT='(
     >  '' Tag filter (P[-rimary]/S[-econdary]/*[any] ; '',
     >  ''now : '',A1,'') : '' )')  KLET
        READ(5,FMT='(A1)',ERR=52) KLET
        IF(KLET.NE.'P' .AND. KLET.NE.'S' .AND. KLET.NE.'*') GOTO 52
        CALL READCA(KKEX,KLET)
      GOTO 98

 6    CONTINUE
C------- Select pass number
        CALL READC2(5)
      GOTO 98

 7    CONTINUE
C------- Select lmnt number
C        LM0 = LM
C 71     WRITE(6,107) 
C 107    FORMAT('   Give # of the optical element to be ray-traced in',
C     >       /,'   (-1 for all)')
C        READ(5,*,ERR=71) LM
C        IF(LM .NE. LM0) OKECH = .FALSE.
        CALL READC4(5)
        CALL READC3(KL1,KL2)
        LM = KL1
      GOTO 98

 8    CONTINUE
C------- Select particle number
        CALL READC6(5)
      GOTO 98

 10   CONTINUE
        CALL TYPTRA
      GOTO 98

 11   CONTINUE
        CALL PLOTE3
      GOTO 98

 12   CONTINUE
        CALL PLOTE1
      GOTO 98

 13   CONTINUE
        IF(OKXAV) THEN
          IF(OKYAV) THEN
            IOP=3
          ELSE
            IOP=1
          ENDIF
        ELSE
          IF(OKYAV) THEN
            IOP=2
          ELSE
            IOP=0
          ENDIF          
        ENDIF
 131    WRITE(6,FMT='('' Give your choice :  y vs. x,  y vs. x-<x>,  '',
     >    ''y-<y> vs. x,  y-<y> vs. x-<x>, (0/1/2/3)'',/,
     >   ''   (now : '',I1,'') : '')') IOP
        READ(5,*,ERR=131) IOP
        IF(IOP.LT.0 .OR. IOP .GT. 3) GOTO 131
        OKXAV = (IOP .EQ. 1) .OR.  (IOP .EQ. 3) 
        OKYAV = (IOP .EQ. 2) .OR.  (IOP .EQ. 3) 
        CALL PLOT4(OKXAV)
        CALL PLOT41(OKYAV)
      GOTO 98

 14   CONTINUE
      WRITE(6,147) 
 147  FORMAT(/,5X,' Possible changes to quantities to be plotted : ', 
     >/,8X,' 1:  horiz. axis  :   X -> a*X^n  + b   ', 
     >/,8X,' 2:  horiz. axis  :   X -> a*(X-b)^n    ', 
     >/,8X,' 3:  horiz. axis  :   X -> a*|X|^n + b   ', 
     >/,8X,' 4:  horiz. axis  :   X -> a*|X-b|^n   ', 
     >/,8X,' 5:  vertic. axis :   Y -> c*Y^m  + d   ', 
     >/,8X,' 6:  vertic. axis :   Y -> c*(Y-d)^m ', 
     >/,8X,' 7:  vertic. axis :   Y -> c*|Y|^m + d   ', 
     >/,8X,' 8:  vertic. axis :   Y -> c*|Y-d|^m   ', 
     >/,8X,'88:  RESET ', 
     >/)
      WRITE(6,148) ' * Option  number : '
 148  FORMAT(A20)
      READ(5,*,ERR=14) IO

      IF(IO.EQ.88) GOTO 1488
      GOTO(141,141,141,141,142,142,142,142) IO
      GOTO 14

 141    WRITE(6,FMT='('' H axis, give your choice for   a, n, b : '', 
     >   /,''    (b=-99 if read from zpop_7.3.14.in)'', 
     >   /,'' present values : '',1P,3G12.4,'' : '')') AX,PX,BX
        READ(5,*,ERR=141)  AX,PX,BX
        WRITE(6,FMT='('' a, n, b : '', 1P,3G12.4,'') : '')') 
     >        AX,PX,BX
        IBXY = 0
        IF(BX.EQ.-99.D0) IBXY = 1
        GOTO 149

 142    WRITE(6,FMT='('' V axis, give your choice for   c, m, d : '', 
     >   /,''    (b=-99 if read from zpop_7.3.14.in)'', 
     >   /,'' present values : '',1P,3G12.4,'' : '')') AY,PY,BY
        READ(5,*,ERR=142)  AY,PY,BY
        WRITE(6,FMT='('' c, m, d : '', 1P,3G12.4,'') : '')') 
     >        AY,PY,BY
        IBXY = 0
        IF(BY.EQ.-99.D0) IBXY = 2
        GOTO 149

 1488 CONTINUE
        AX = 1.D0
        PX = 1.D0
        BX = 0.D0
        AY = 1.D0
        PY = 1.D0
        BY = 0.D0
        IBXY = 0
        CALL PLOT51(AX,PX,BX,AY,PY,BY)
        CALL BIN41W(AX,PX,BX,AY,PY,BY)
      GOTO 98

 149  CONTINUE
        IF(IBXY.NE.0) THEN
          IF (IDLUNI(IUN)) THEN
            OPEN(UNIT=IUN,FILE='zpop_7.3.14.in')
          ELSE
            WRITE(6,*) 
     >      ' *** SBR PLTOPT, problem : could not open  zpop_7.3.14.in'
          ENDIF
        ENDIF
        CALL PLOT5(AX,PX,BX,AY,PY,BY,IO,IBXY,IUN)
        CALL BIN4W(AX,PX,BX,AY,PY,BY,IO,IBXY,IUN)
      GOTO 98

 15   CONTINUE
        RETURN
C------- Plot X(ipass)-X(ipass-1), Y(ipass)-Y(ipass-1)
        IF(OKX12) THEN
          IF(OKY12) THEN
            IOP=3
          ELSE
            IOP=1
          ENDIF
        ELSE
          IF(OKY12) THEN
            IOP=2
          ELSE
            IOP=0
          ENDIF          
        ENDIF
 151    WRITE(6,FMT='('' Give your choice :  y vs. x,  y vs. dx,  '',
     >    ''dy vs. x,  dy vs. dx, (0/1/2/3)'',/,
     >   ''   (now : '',I1,'') : '')') IOP
        READ(5,*,ERR=151) IOP
        IF(IOP.LT.0 .OR. IOP .GT. 3) GOTO 151
        OKX12 = (IOP .EQ. 1) .OR.  (IOP .EQ. 3) 
        OKY12 = (IOP .EQ. 2) .OR.  (IOP .EQ. 3) 
      GOTO 98

 26   CONTINUE
C------- Path extrapolation 
        CALL PLOTE2
      GOTO 98

 88   RETURN
 98   RETURN
      END
