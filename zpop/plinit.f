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
C  Brookhaven National Laborator
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE PLINIT(TXTK,
     >                       NL,OKOPN,CHANGE,FILE2O)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) TXTK, FILE2O
      LOGICAL OKOPN,CHANGE

      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN
      COMMON/VEN/Y1,S1,D,DQ,DS,DM,DSS,DPU,TETA1,
     >           SMIN,SMAX,YMIN,YMAX,SB0,YB0
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      CHARACTER(*) PART, GNUFIL
      DIMENSION OX(*), WF(*), WO(*), FNR(*)
      LOGICAL STORE

      CHARACTER(80) TITL, FILE3O
      CHARACTER(200) TXT
      CHARACTER(5) STOR
      LOGICAL EXS, OPN, IDLUNI 
      LOGICAL EXS2, OPN2
      INTEGER DEBSTR
      LOGICAL EMPTY
      SAVE LUO
      LOGICAL TYLAB

      parameter (mst=9)
      character(20) stra(mst), tmod
      logical strcon
      
      DATA IFIVE /5/
      
      JUN = 6
      CALL IMPV2(JUN)

      ENTRY SRINIT(GNUFIL,TXTK,
     >                         PART,R0,Q,AM,OX,WF,WO,FNR,STORE)
      kxo = kx ; kyo = ky
      xmio = xmi; xmao = xma; ymio = ymi; ymao = yma;
      liso = lis ; nbo = nb
      
      ENTRY PLINI1
c      ENTRY PLINI1(TXTK,
c     >     kxo,kyo,tmod,xco,xpco,alf,bet,
c     >     xmio,xmao,ymio,ymao,LISo, NBo,mark,
c     >      FILE2O, kp1,kp2,kp3,kt1,kt2,FILE3O )

c      kx = kxo ; ky = kyo
c      xmi = xmio; xma = xmao; ymi = ymio; yma = ymao;
c      lis = liso ; nb = nbo
      
      INQUIRE(FILE='zpop.dat',EXIST=EXS, ERR=98,IOSTAT=IOS,
     >  OPENED=OPN, NUMBER=LN)
      IF(EXS) THEN
        IF(OPN) THEN
          REWIND(LN)
        ELSE
          IF (IDLUNI(LN)) THEN
            OPEN(UNIT=LN,FILE='zpop.dat',ERR=99,IOSTAT=IOS)
            IF(IOS.NE.0) GOTO 99
          ELSE
            WRITE(6,*) ' *** SBR PLINIT ; Error idle UNIT '
            RETURN 
          ENDIF
        ENDIF

  5     CONTINUE
          READ(LN,FMT='(A)',ERR=92,END=95) TXT
          IDB=DEBSTR(TXT)
          IF(TXT(IDB:IDB+LEN(TXTK)-1).NE.TXTK) GOTO 5

        READ(LN,FMT='(A)',ERR=95,END=95) TITL

        IF(TXTK.EQ.'////MENU08/16////') THEN
C---------- Synchrotron radiation menu
          STOR = ' '
          READ(LN,FMT='(A)',ERR=12,END=10) PART
          READ(LN,*,ERR=12,END=10) R0, Q, AM
          READ(LN,*,ERR=12,END=10) (OX(I), I=1, 3)
          READ(LN,*,ERR=12,END=10) (WF(I), I=1, 3)
          READ(LN,*,ERR=12,END=10) (WO(I), I=1, 6)
          READ(LN,*,ERR=12,END=10) (FNR(I), I=1, 3)
          READ(LN,FMT='(A)',ERR=12,END=10) STOR

 10       CONTINUE
          STORE = (STOR .EQ. 'store' .OR. STOR .EQ. 'STORE')
          CALL SRINM1(STORE)
          IF(STORE) THEN 
            WRITE(*,*) ' Pgm plinit. STORE option found  ',STOR,
     >      '  in zpop.dat'
            CALL OPNGNU(GNUFIL, 
     >                                 STORE,LUO)
          ENDIF
          WRITE(6,*) ' Parameters for SR calculation read from ',
     >                  'the data file "zpop.dat", '
          WRITE(6,FMT='(2A)') '    - this data list is entitled : ',TITL
          GOTO 11

 12       WRITE(6,*) 
          WRITE(6,*) ' Premature end of read in data file "zpop.dat" ! '
          WRITE(6,*) '         Probably missing data... '

 11       CONTINUE

        ELSEIF(TXTK.EQ.'////MENU07////') THEN
C---------- Main menu
          KX0 = KX
          KY0 = KY
C--------- Lab coordinate types
          LINE = 0
C          READ(LN,*,ERR=57,END=57) KX,KY
          READ(LN,FMT='(A)',ERR=92,END=95) TXT
          if(strcon(txt,'!',
     >                      is)) txt = txt(1:is-1)
          call strget(txt,9,
     >                      nst,stra)
          if(nst .gt. mst) nst=mst
          
c          write(*,*) ' Pgm plinit txt ',txt
c            write(*,*) ' Pgm plinit stra ',nst,(stra(i),' ',i=1,nst)
c            read(*,*) 
            
          if(nst .ge. 2) then
            read(stra(1),*,err=92,end=92) KX
            read(stra(2),*,err=92,end=92) KY
            if(nst .ge. 3) then
              tmod = stra(3)
              if(tmod .eq. 'N') then
                if(nst .ge. 9) then
                  read(stra(4),*,err=92,end=92) xco
                  read(stra(5),*,err=92,end=92) xpco 
                  read(stra(6),*,err=92,end=92) alf
                  read(stra(7),*,err=92,end=92) bet               
                  read(stra(8),*,err=92,end=92) eta
                  read(stra(9),*,err=92,end=92) etap               
                else
                  write(*,*) ' Pgm plinit. Error input alfa, beta.'
                  goto 92
                endif
              elseif(tmod .eq. 'normD') then

              else
                tmod = ' '
              endif
            endif
          else
            write(*,*) ' Pgm plinit. Error input kx, ky.'
            goto 92
          endif

c            write(*,*) ' stra(i) ',nst,(stra(i),' ',i=1,nst)
c            read(*,*)
            
          WRITE(*,*) ' Pgm plinit. KX, KY : ',KX,KY
          if(nst .ge. 3) WRITE(*,*)
     >         ' Pgm plinit. nst, alfa, beta, eta, etap : ',
     >          nst,xco,xpco,alf,bet,eta,etap
          call plot20(tmod,xco,xpco,alf,bet,eta,etap)
          LINE = LINE+1
          OKVAR = .TRUE.
          IF( TYLAB(KX,KY) ) THEN
              IF(KX.NE.48) THEN
                KY=KX
                KX=48
              ENDIF
c          ELSEIF(KX.EQ.68 .AND. KY.EQ.62) THEN
c            MOD = 0
c            RM = 3.482590E2 
c            CM2M = 0.01
c            CALL READCC(MOD,RM*CM2M,RM*CM2M)
          ENDIF
C--------- Axis range
          XMI0 = XMI
          XMA0 = XMA
          YMI0 = YMI
          YMA0 = YMA
          READ(LN,*,ERR=57,END=57) XMI,XMA,YMI,YMA
          WRITE(*,*) ' Pgm plinit. XMI,XMA,YMI,YMA : ',
     >    XMI,XMA,YMI,YMA
          LINE = LINE+1
          IF(XMI.LT. XMA .AND. YMI.LT. YMA) THEN
            CALL TRAXES(XMI,XMA,YMI,YMA,2)
            OKECH=.TRUE.
          ELSE
            XMI = XMI0
            XMA = XMA0
            YMI = YMI0
            YMA = YMA0
            OKECH = .FALSE.
          ENDIF
C--------- LIST option and histogram bins number
          READ(LN,*,ERR=57,END=57) LIS, NB
          WRITE(*,*) ' Pgm plinit. LIS, NB : ',LIS, NB
          LINE = LINE+1
C--------- Marker/line type
          READ(LN,*,ERR=57,END=57) MARK
          WRITE(*,*) ' Pgm plinit. MARK : ', MARK
          LINE = LINE+1
          CALL LINTYW(MARK)
C--------- Trajectories data file (default is zgoubi.plt)
          FILE2O = ' '
          READ(LN,*,ERR=57,END=57) FILE2O
          WRITE(*,*) ' Pgm plinit. FILE2O : ',FILE2O
          LINE = LINE+1
          READ(LN,*,ERR=57,END=57) KP1, KP2, KP3
          IF(KP1.EQ.-1) THEN
            WRITE(*,*) ' Pgm plinit. Old style KP1-3, '
     >      //' converted to new style '  
            KP1 = 1
            KP3 = KP2
            KP2 = 999999
          ENDIF
          WRITE(*,*) ' KP1-3 : ',KP1, KP2, KP3
          LINE = LINE+1
          READ(LN,*,ERR=57,END=57) KT1, KT2
          WRITE(*,*) KT1, KT2
          LINE = LINE+1
C--------- Optical structure data file (default is zgoubi.dat)
          FILE3O = ' '
          READ(LN,*,ERR=57,END=57) FILE3O  !!!!!, LMNT1, LMNT2
          WRITE(*,*) FILE3O  !!!!!, LMNT1, LMNT2
          LINE = LINE+1
C--------- Options, menu 7.3.14
          READ(LN,*,ERR=57,END=57) AX,PX,BX
          WRITE(*,*) AX,PX,BX
          IBXY = 0
          IF(BX.EQ.-99.D0) IBXY = 1
          IF(IBXY.NE.0) THEN
            INQUIRE(FILE='zpop_7.3.14.in',EXIST=EXS2, ERR=98,IOSTAT=IOS,
     >      OPENED=OPN2, NUMBER=LN2)
            IF(EXS2) THEN
              IF(OPN2) THEN
                REWIND(LN2)
              ELSE
                IF (IDLUNI(IUN)) THEN
                  OPEN(UNIT=IUN,FILE='zpop_7.3.14.in')
                  OPN2 = .TRUE.
                ELSE
                  WRITE(6,*) 
     >        '*** SBR PLINIT, problem : could not open  zpop_7.3.14.in'
                ENDIF
              ENDIF
            ELSE
              WRITE(6,*) 
     >        '*** SBR PLINIT : no such file to open, zpop_7.3.14.in'
            ENDIF
          ENDIF
          CALL PLOT5(AX,PX,BX,DX,AY,PY,BY,DY,IO,IBXY,IUN)
          CALL BIN4W(AX,PX,BX,DX,AY,PY,BY,DY,IO)
          LINE = LINE+1
          READ(LN,*,ERR=57,END=57) AY,PY,BY
          WRITE(*,*) AY,PY,BY
          IBXY = 0
          IF(BY.EQ.-99.D0) IBXY = 2
          IF(IBXY.NE.0) THEN
            INQUIRE(FILE='zpop_7.3.14.in',EXIST=EXS2, ERR=98,IOSTAT=IOS,
     >      OPENED=OPN2, NUMBER=LN2)
            IF(EXS2) THEN
              IF(OPN2) THEN
                REWIND(LN2)
              ELSE
                IF (IDLUNI(IUN)) THEN
                  OPEN(UNIT=IUN,FILE='zpop_7.3.14.in')
                  OPN2 = .TRUE.
                ELSE
                  WRITE(6,*) 
     >        '*** SBR PLINIT, problem : could not open  zpop_7.3.14.in'
                ENDIF
              ENDIF
            ELSE
              WRITE(6,*) 
     >        '*** SBR PLINIT : no such file to open, zpop_7.3.14.in'
            ENDIF
          ENDIF
          CALL PLOT5(AX,PX,BX,DX,AY,PY,BY,DY,IO,IBXY,IUN)
          CALL BIN4W(AX,PX,BX,DX,AY,PY,BY,DY,IO)
          LINE = LINE+1
C--------- Options, menu 7.7(LabCoord).5. Initial coordinates of synoptic
          READ(LN,*,ERR=57,END=57) SB0, YB0, TETA1   ! (m, m, rad)
          WRITE(*,*) SB0, YB0, TETA1
          CALL INSY4(SB0,YB0,TETA1)
          LINE = LINE+1
          
 57       WRITE(6,FMT='(I2,2A)') LINE,' lines of initialisation ',  
     >      ' parameters have been read from the data file "zpop.dat" '
          WRITE(6,FMT='(2A)') '    - this data list is entitled : ',TITL
          IF(.NOT.EMPTY(FILE2O)) THEN
            IFIVE=5
            CALL OPNMN(IFIVE,
     >                       NL,OKOPN,CHANGE,FILE2O)
            WRITE(6,FMT='(A,/)')' for reading coordinates to be plotted'
          ENDIF
          IF(.NOT.EMPTY(FILE3O)) CALL INSY2(FILE3O)
          CALL READC2B(KP1,KP2,KP3)
          CALL READC6B(KT1,KT2)
        ENDIF

        CLOSE(LN)
        RETURN        

      ELSE
        WRITE(6,*) 
        WRITE(6,*) ' Data file "zpop.dat" does not exist. '
        GOTO 97
      ENDIF

 92   WRITE(6,*) ' SBR PLINIT : Error READ in zpop.dat at LINE=',line
      GOTO 97
 95   WRITE(6,*) ' SBR PLINIT : End of file "zpop.dat" reached, ', 
     >     ' could not find initialisation data.'
      GOTO 97
 99   WRITE(6,*) 
     >  ' SBR PLINIT : Warning, Cannot open data file "zpop.dat" ; '
      GOTO 97
 98   WRITE(6,*) 
     >  ' SBR PLINIT : Warning, INQUIRE(FILE=''zpop.dat'') failed ;'
      WRITE(6,*)
     >  '           Check existence of zpop.dat.'
 97   WRITE(6,*) '        Default initialisation values will be taken '

      IF(TXTK.EQ.'////MENU08/16////') THEN
        WRITE(6,*) '   for power normalisation, particle type,', 
     >                                  ' observer position, etc.'
        PART = 'E'
        IF    (PART .EQ. 'P') THEN
          AM = .93827231D3
          R0 = 1.53469852D-18
          Q = 1.60217733D-19
        ELSE       !!!IF(PART .EQ. 'E') THEN
          AM =  .51099906D0
          R0 =  2.81794092D-15
          Q = 1.60217733D-19
        ENDIF      
C------- Observer position, x, y, z (m)
        OX(1) = 1.D5
        OX(2) = 0.D0
        OX(3) = 0.D0
C------- Frequency range for spectrum calculation
C        Nu1, Nu2 (keV), number of samples
        WF(1) = 1.D-9
        WF(2) = 1.D0
        WF(3) = 999.D0
C------- Window width : dummy, +/- WY, +/- Wz,
C                     dummy, fragmentation NY and NZ for integration
        WO(1) = 0.D0
        WO(2) = 10.D0
        WO(3) = 10.D0
        WO(4) = 0.D0
        WO(5) = 7.D0
        WO(6) = 7.D0
C------- The spectrum and integrals are normalized by FNR(1)/FNR(2)*FNR(3)
C      Normally FNR(1)/FNR(2) = I(A)/q(C) and FNR(3) = 1,
C        in order to get the power (Watt)  emitted by a circulating beam, 
C        from the energy radiated by 1 particle, 1 pass.
        FNR(1) = 1.D0
        FNR(2) = 1.D0
        FNR(3) = 1.D0
      ELSEIF(TXTK.EQ.'////MENU07////') THEN
C------- Initialisation data are in SBR BLOCK
        WRITE(6,*) '   for variables to plot, line type, pass #, etc.'
      ENDIF
      RETURN
      END
