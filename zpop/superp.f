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
      SUBROUTINE SUPERP(OKECH,NLOG,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKECH

      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB            

      LOGICAL EMPTY, IDLUNI, STRCON
      SAVE LN
      CHARACTER TXTX*15, TXTY*15, TITL*70, FNAMI*120, FNAM*120, TXT*80, 
     >            TXT80*80,TXT4*4
      SAVE FNAM
      PARAMETER (MTB=20)
      DIMENSION TAB(MTB)
      PARAMETER (PTXT0=14.D0)
      SAVE PTXT
      CHARACTER * 9   DMY
      CHARACTER * 1 CTR, CTR0
      CHARACTER * 1 REP
      SAVE MOD
C      SAVE XN,YN
      INTEGER DEBSTR, FINSTR

      SAVE AX,PX,BX,AY,PY,BY

      DATA FNAM / 'fort.88' /
      DATA IX, JY / 1, 2 /
      DATA PTXT / PTXT0 /
      DATA CTR0 / 'N' /     
C      DATA XN, YN / 1.D0, 1.D0 /
      DATA MOD / 2 /
      DATA DX / 0.D0 /
      DATA AX,PX,BX,AY,PY,BY / 1.D0, 1.D0, 0.D0, 1.D0, 1.D0, 0.D0 /

      DATA I0,ZERO / 0, 0.D0/

C      LTYP = LTYP + 1
C      IF(LTYP .EQ. 6) LTYP = 1
C      CALL LINTYW(LTYP)

      CALL LINTYR(LTYP0)  
      LTYP = LTYP0

      WRITE(6,*) 
     > '  ** Your input data file is expected to have no header **'
      WRITE(6,*) '   It  is expected to contain data in columns'
      WRITE(6,*) '     - MAXIMUM ALLOWED IS ',MTB,' DATA COLUMNS -'
      WRITE(6,*) '   You can plot any of these, against any other'
      WRITE(6,*) 
     >  '   Blank lines are considered as separators between curves'
      WRITE(6,*)
 20   WRITE(6,*) ' Give the name of the input data file'
      WRITE(6,*) '           - default will be ',FNAM
      READ(5,FMT='(A)',ERR=20) FNAMI
      IF( .NOT. EMPTY(FNAMI) ) THEN
        FNAMI = FNAMI(DEBSTR(FNAMI):FINSTR(FNAMI))
        FNAM = FNAMI
      ENDIF
      IF(LN.GT.0) CLOSE(LN)
      IF (IDLUNI(LN)) THEN
        WRITE(6,*) ' SBR SUPERP, Logical unit ',ln,' is free...'
        WRITE(6,*) ' SBR SUPERP, Try to open ',FNAM
        OPEN(UNIT=LN,FILE=FNAM,STATUS='OLD',ERR=98)
        WRITE(6,*) ' Opened input data file ',FNAM
      ELSE
        WRITE(6,*) ' *** SBR SUPERP ; Error idle UNIT '
        RETURN 1
      ENDIF

 21   WRITE(6,*) ' Plot column numbers (x & y) - default is ',IX,JY,':'
      READ(5,FMT='(A)',ERR=21) TXT
      READ(TXT,*,ERR=210,END=210) IXI, JYI
 210  IF(IXI.GE.1 .AND. IXI .LE. MTB) IX=IXI
      IF(JYI.GE.1 .AND. JYI .LE. MTB) JY=JYI
      WRITE(6,*) ' Will  plot  col.# ',JY,'  vs.  col.# ',IX
      WRITE(6,*) 

      WRITE(6,*) ' x and y norm.  (y/n) : '
      READ(5,FMT='(A)',ERR=26) REP
      IF( EMPTY(REP) ) REP = 'n'
      IF(REP .EQ. 'y') REP = 'Y'
      IF(REP .NE. 'Y') REP = 'N'
      IF(REP.EQ.'Y') THEN
 141     WRITE(6,FMT='('' Give your choice for a*x^p+b,   a, p, b : '', 
     >   /,''    (now : '',1P,3G12.4,'') : '')') AX,PX,BX
          READ(5,*,ERR=141)  AX,PX,BX
C          PX = NINT(XN)
          WRITE(6,FMT='('' a, n, b : '', 1P,3G12.4,'') : '')') 
     >          AX,PX,BX
 142      WRITE(6,FMT='('' Give your choice for c*y^q+d,   c, q, d : '', 
     >     /,''    (now : '',1P,3G12.4,'') : '')') AY,PY,BY
          READ(5,*,ERR=142)  AY,PY,BY
C          PY = NINT(YN)
          WRITE(6,FMT='('' c, m, d : '', 1P,3G12.4,'') : '')') 
     >        AY,PY,BY
      ENDIF

 27   WRITE(6,*) ' Lin/lin axis, log/log axis (2/32) (now ',MOD,') :'
      READ(5,FMT='(A)',ERR=27) TXT
      READ(TXT,*,ERR=270,END=270) MODR
      IF(MODR .NE. 1 .AND. MODR .NE. 31
     >   .AND. MODR .NE. 2 .AND. MODR .NE. 32 ) GOTO 27
      MOD = MODR      
 270  CONTINUE

 26   WRITE(6,*) ' Center x coord on graph (n/y) :'
      READ(5,FMT='(A)',ERR=26) CTR
      IF( EMPTY(CTR) ) CTR = CTR0
      CTR0=CTR
      IF(CTR .EQ. 'y') CTR = 'Y'
      IF(CTR .NE. 'Y') CTR = 'N'

 22   WRITE(6,*) ' Give X-axis text - Max 15 char. Default is blank :'
      READ(5,FMT='(A)',ERR=22) TXTX
      IF( EMPTY(TXTX) ) TXTX = ' '
 23   WRITE(6,*) ' Give Y-axis text - Max 15 char. Default is blank :'
      READ(5,FMT='(A)',ERR=23) TXTY
      IF( EMPTY(TXTY) ) TXTY = ' '
 24   WRITE(6,*)' Give figure subtitle - Max 80 char. Default is blank:'
      READ(5,FMT='(A)',ERR=24) TITL
      IF( EMPTY(TITL) ) TITL = ' '



C       OKECH=.false.
      MXTAB = IX
      IF(JY .GT. IX) MXTAB = JY

      IF( .NOT. OKECH ) THEN
        XMI = 1.D10
        XMA =-1.D10
        YMI = 1.D10
        YMA =-1.D10
        NR = 1

C        READ(LN,*,END=10,ERR=3) (TAB(I), I=1, MXTAB)
 11     READ(LN,FMT='(A80)',END=10,ERR=3) TXT80
        IF    (STRCON(TXT80,'%',
     >                       IS) ) THEN
          GOTO 11
        ELSEIF(STRCON(TXT80,'&',
     >                       IS) ) THEN
          GOTO 11
        ELSEIF(EMPTY(TXT80)) THEN
          GOTO 11
        ELSEIF(STRCON(TXT80,'LTYP',
     >                             IS) ) THEN
          GOTO 11
        ELSE
          BACKSPACE(LN)
          READ(LN,*,ERR=3) (TAB(I), I=1, MXTAB)
        ENDIF

c          write(*,*) (TAB(I), I=1, MXTAB), nr
        XX = TAB(IX) 
CCCCC-------------- Rustine pour plot IR1 LHC
CCCCC           if(x.lt. 500) x = x+2.6659E+04
        YY = TAB(JY) 

        X = AX*XX**PX + BX
        Y = AY*YY**PY + BY

        CALL MINMAX(X,Y,XMI,XMA,YMI,YMA)

 1      CONTINUE
          NR = NR+1
C          READ(LN,*,END=2,ERR=3) (TAB(I), I=1, MXTAB)
          READ(LN,FMT='(A80)',END=2,ERR=3) TXT80
          IF    (STRCON(TXT80,'%',
     >                         IS) ) THEN
            GOTO 1
          ELSEIF(STRCON(TXT80,'&',
     >                       IS) ) THEN
            GOTO 1
          ELSEIF(EMPTY(TXT80)) THEN
            GOTO 1
          ELSEIF(STRCON(TXT80,'LTYP',
     >                               IS) ) THEN
            GOTO 1
          ELSE
            BACKSPACE(LN)
            READ(LN,*,ERR=3) (TAB(I), I=1, MXTAB)
          ENDIF

c        write(*,*) (TAB(I), I=1, MXTAB), nr

          XX = TAB(IX)
CCCCCC-------------- Rustine pour plot IR1 LHC
CCCCC           if(x.lt. 500) x = x+2.6659E+04
          YY = TAB(JY)

          X = AX*XX**PX + BX
          Y = AY*YY**PY + BY

          CALL MINMAX(X,Y,XMI,XMA,YMI,YMA)
        GOTO 1

 3      CONTINUE
        CALL FBGTXT
        WRITE(6,*) ' *** SBR SUPERP ; Error on READ in ',FNAM
        WRITE(6,*) '                  at event # ', NR
        WRITE(6,*) '          during axis min-max calculation'
 2      CONTINUE

        IF(XMI.LT. XMA .AND. YMI.LT. YMA) THEN
          DX = 0.D0
          IF(CTR .EQ. 'Y') THEN
            DX = (XMA+XMI)/2.D0
            XMA = XMA-DX
            XMI = -XMA
          ENDIF
C          DDX=ABS(XMA-XMI)*0.05D0
C          DDY=ABS(YMA-YMI)*0.05D0
          CALL TRAXES(XMI,XMA,YMI,YMA,MOD) 
C          CALL TRAXES(XMI-DDX,XMA+DDX,YMI-DDY,YMA+DDY,MOD)
          OKECH = .TRUE.
        ELSE
          WRITE(6,*) ' Xmin-max, Ymin-max :',XMI,XMA,YMI,YMA
          GOTO 96
        ENDIF
      ENDIF
C----- endif okech

C----- NOW PLOT !
      REWIND(LN)

      CALL LINTYP(-1)

      XMIL = 1.D10
      XMAL =-1.D10
      YMIL = 1.D10
      YMAL =-1.D10
      NR = 1
      NPT = 0

C      READ(LN,*,END=10,ERR=33) (TAB(I), I=1, MXTAB)
 12   READ(LN,FMT='(A80)',END=10,ERR=33) TXT80
      IF    (STRCON(TXT80,'%',
     >                     IS) ) THEN
        GOTO 12
      ELSEIF(STRCON(TXT80,'&',
     >                       IS) ) THEN
        GOTO 12
      ELSEIF(EMPTY(TXT80)) THEN
        GOTO 12
      ELSEIF(STRCON(TXT80,'LTYP',
     >                           IS) ) THEN
        txt80 = txt80(debstr(txt80):80)
        READ(TXT80,FMT='(a,i2)') txt4,LTYP
        CALL LINTYP(LTYP)
        MODV=4
        GOTO 12
      ELSE
        BACKSPACE(LN)
        READ(LN,*,ERR=33) (TAB(I), I=1, MXTAB)
      ENDIF

c          write(*,*) (TAB(I), I=1, MXTAB), nr
      XX = TAB(IX)
CCCCCC-------------- Rustine pour plot IR1 LHC
CCCCCC           if(x.lt. 500) x = x+2.6659E+04
      YY = TAB(JY)

          X = AX*XX**PX + BX - DX
          Y = AY*YY**PY + BY

      CALL MINMAX(X,Y,XMIL,XMAL,YMIL,YMAL)
      AVX = X
      AVY = Y
      SQX = X*X
      SQY = Y*Y

      CALL VECTPL(X,Y,4)
      CALL VECTPL(X,Y,2)
      NPT = 1
C      X2 = X
      MODV=2

 4    CONTINUE
        NR = NR+1
        READ(LN,FMT='(A80)',END=10,ERR=33) TXT80
        IF    (STRCON(TXT80,'%',
     >                      IS) ) THEN
          MODV=4
          GOTO 4
        ELSEIF(STRCON(TXT80,'&',
     >                       IS) ) THEN
          MODV=4
          GOTO 4
        ELSEIF(EMPTY(TXT80)) THEN
          LTYP = LTYP+1
          IF(LTYP .EQ. 6) LTYP = 1
          CALL LINTYP(LTYP)
          MODV=4
          GOTO 4
        ELSEIF(STRCON(TXT80,'LTYP',
     >                             IS) ) THEN
          txt80 = txt80(debstr(txt80):80)
          READ(TXT80,FMT='(a,i2)') txt4,LTYP
          CALL LINTYP(LTYP)
          MODV=4
          GOTO 4
        ELSE
          BACKSPACE(LN)
          READ(LN,*,ERR=33) (TAB(I), I=1, MXTAB)
        ENDIF

c          write(*,*) (TAB(I), I=1, MXTAB), nr
        XX = TAB(IX)
CCCCCCC-------------- Rustine pour plot IR1 LHC
CCCCCC           if(x.lt. 500) x = x+2.6659E+04
        YY = TAB(JY)

          X = AX*XX**PX + BX - DX
          Y = AY*YY**PY + BY

        CALL MINMAX(X,Y,XMIL,XMAL,YMIL,YMAL)
        AVX = AVX + X
        AVY = AVY + Y
        SQX = SQX + X*X
        SQY = SQY + Y*Y

        CALL VECTPL(X,Y,MODV)
        IF(MODV.EQ.4) MODV=2

        NPT=NPT+1
        IF(LIS .EQ. 2)CALL IMPV(NLOG,NPT,X,Y,ZERO,ZERO,I0)

      GOTO 4

 33   CONTINUE
      CALL FBGTXT
      WRITE(6,*) ' *** SBR SUPERP ; end plot on error READ in ',FNAM
      WRITE(6,*) '                  at event # ', NR
      NR = NR - 1

 10   CONTINUE
C      CALL VECTPL(X,Y,4)
      CALL VECTPL(X,Y,2)
C      NPT=NPT+1
      CLOSE(LN)
      CALL FBGTXT

      AVX = AVX/NPT
      AVY = AVY/NPT
      SQX = SQX / NPT
      SQY = SQY / NPT
      SIGX = SQRT(SQX - AVX*AVX)
      SIGY = SQRT(SQY - AVY*AVY)

      WRITE(6,100) TXTY, TXTX
 100  FORMAT(A,'vs. ',A) 
      WRITE(6,107) TITL
 107  FORMAT(5X,A)
      WRITE(6,101) AVX, SIGX, XMIL, XMAL
 101  FORMAT(5X,'<X>, Sig_X, X_min,_max : ',1P,4G12.4)
      WRITE(6,102) AVY, SIGY, YMIL, YMAL
 102  FORMAT(5X,'<Y>, Sig_Y, Y_min,_max : ',1P,4G12.4)
      WRITE(6,FMT='(/,I6,'' points have been plotted'')') NPT
      WRITE(6,*)
      WRITE(6,*) '<X>, Sig_X, <Y>, Sig_Y, X_min, _max, Y_min, _max : '
      WRITE(6,FMT='(1P,8G12.4,/)') AVX,SIGX,AVY,SIGY,XMIL,XMAL,YMIL,YMAL

      CALL TXTFBG
      CALL DATE2(DMY)
      WRITE(TXT,FMT='(A9)') DMY
      CALL DEFCAR(1,0,0)
      CALL TRTXT(38.D0,245.D0,TXT,50,0)
      CALL DFKSIZ
      WRITE(TXT,100) TXTY, TXTX
      CALL TRTXT(120.D0,245.D0,TXT,80,0)
      WRITE(TXT,FMT='(A54)') '* '//TITL(1:50)//' *'
      CALL TRTXT(60.D0,28.D0,TXT,50,0)
      WRITE(TXT,101) AVX, SIGX, XMIL, XMAL
      CALL TRTXT(10.D0,PTXT,TXT,80,0)
      PTXT = PTXT - 7.D0
      WRITE(TXT,102) AVY, SIGY, YMIL, YMAL
      CALL TRTXT(10.D0,PTXT,TXT,80,0)
      PTXT = PTXT - 7.D0
      WRITE(TXT,FMT='(I6,'' points plotted'')') NPT
      CALL TRTXT(20.D0,PTXT,TXT,80,0)
      IF(PTXT .LT. 0.D0) PTXT = PTXT0

      CALL FBGTXT
      CLOSE(LN)
      RETURN

 96   CONTINUE
      CALL FBGTXT
      WRITE(6,*) ' *** SBR SUPERP ; Error on axis Min-Max '
      CLOSE(LN)
      RETURN 1
 98   CONTINUE
      WRITE(6,*) ' *** SBR SUPERP ; Error on OPEN ',FNAM
      CLOSE(LN)
      RETURN 1

      END
