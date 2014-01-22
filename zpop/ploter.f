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
      SUBROUTINE PLOTER(NLOG,NL,KPS,NPTS,NPTR,OKREW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------------------
C     TRACE LES VARIABLES KY V.S. KX 
C     ------------------------------
      LOGICAL OKREW
      PARAMETER (MXB=1000)
      COMMON/B/BINS(MXB+2)
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN
      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*7, KPOL(2)*9, KDIM(MXVAR)*7
      COMMON/INPVR/ KVAR, KPOL, KDIM
      COMMON/LUN/NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB
      COMMON/REBELO/ NRBLT,IPASS,KWRI,NNDES,STDVM
      
      CHARACTER LET0
      CHARACTER LET
      DIMENSION YZXB(MXVAR),NDX(5)
      SAVE YZXB
      LOGICAL INRANG

      INCLUDE 'MAXNTR.H'
      INCLUDE 'MAXTRA.H'

      DIMENSION XYLAST(2,MXT)
      LOGICAL CNECT
      CHARACTER*3 TXT
      DIMENSION OKT(MXT)
      LOGICAL OKTAG
      SAVE TAGX, TAGY, OKTAG
      SAVE PPA

      LOGICAL OKS0, OKS
      SAVE OKS0, OKS, SMAX
      LOGICAL OKXAV, OKYAV, OKXAVI, OKYAVI, OKXAVO, OKYAVO
      SAVE OKXAV, XMOY, OKYAV, YMOY
      LOGICAL OKX12, OKY12, OKX12I, OKY12I, OKX12O, OKY12O
      SAVE OKX12, OKY12

      INCLUDE "MXSTEP.H"
      DIMENSION BTAB(MXSTEP,2)
      SAVE AX,PX,BX,DX,AY,PY,BY,DY, BTAB

      LOGICAL  OKDX, OKABSX, OKABDX, OKDY, OKABSY, OKABDY
      SAVE     OKDX, OKABSX, OKABDX, OKDY, OKABSY, OKABDY

      LOGICAL BINM

C----- For possible storing of coordinates, if ellipse plot
      LOGICAL KLIPS, KLIPSI                   
      SAVE KLIPS

      LOGICAL IDLUNI
      LOGICAL READBX, READBY

      DATA CNECT / .TRUE. /

C------ Entering negative value for X  will cancel tagging
      DATA TAGX, TAGY, OKTAG / -0.7D0, 1.D-2, .FALSE. /

      DATA     OKDX, OKABSX, OKABDX, OKDY, OKABSY, OKABDY /6*.FALSE./

      DATA KLIPS / .FALSE. /

      DATA OKS / .FALSE./
      DATA PPA / 0.D0 /
      DATA AX,PX,BX,AY,PY,BY / 1.D0, 1.D0, 0.D0, 1.D0, 1.D0, 0.D0 /

      DATA READBX, READBY / .FALSE., .FALSE. /

C      IUNIT = NLOG
C      WRITE(IUNIT,*) '%% <X>, Sig_X, <Y>, Sig_Y, Correlation, IPASS'
C      CALL FLUSH2(IUNIT,.FALSE.)      

      IF(LIS .EQ. 2) THEN 
        IF (IDLUNI(JUN)) THEN
          OPEN(UNIT=JUN,FILE='zpop.log_IMPV')
          CALL IMPV2(JUN)
        ELSE
          WRITE(6,*) '    '
          WRITE(6,*) ' *** sbr PLOTER, Problem : No idle unit number ! '
          WRITE(6,*) '     Forced to 99 '
          WRITE(6,*) '    '
          JUN = 99
          CALL IMPV2(99)
        ENDIF
      ENDIF

      IF(OKREW) CALL REWIN2(NL,*99)

      DO 1 I = 1, MXT
 1      OKT(I) = 0

      NOC = 0
      NRD=0
      D0=0.D0
      ITAB = 1

      LET0 = ' '
      IT0 = 0
      IPASS0 = 0
      CALL LINTYP(-1)
 
      XM = 0.D0
      YM = 0.D0
      X2M = 0.D0
      Y2M = 0.D0
      XYM = 0.D0

      IF(KLIPS) CALL RAZ(BINS,NB)

C----- BOUCLE SUR READ FICHIER NL 
 44   CONTINUE                 
      NRD=NRD+1

C----- Read next coordinate 
      CALL READCO(NL,
     >               KART,LET,YZXB,NDX,*10,*19)
      IPASS=NINT(YZXB(39))
      IF(IPASS .NE. IPASS0) NPT = 0

C----- File type zgoubi.plt
      IF    (NL .EQ. NPLT) THEN

C-------  Particle changes
        IF(LET .NE. LET0 .OR. NDX(2) .NE. IT0 ) THEN

C--------- If this is not the first particle...
C          IF(LET0 .NE. ' ') THEN
          IF(IT0 .NE. 0) THEN

C----------- Extrapolate the plot linearly over an additive path length PPA
C            at end of field maps
            IF(KX.EQ. 8 .AND. KY.EQ. 2..AND.PPA.NE. 0.D0)
     >        CALL PPAR(PPA,X,Y,T0,KART)
          ENDIF

          NPT=0
          itab = 1
          LET0=LET
          IT0 = NDX(2)
           
        ENDIF
        T0 = YZXB(3)

C----- File type zgoubi.fai
      ELSEIF(NL .EQ. NFAI) THEN   

        IF(YZXB(11) .NE. D0) THEN
C---------- Detects a change of dp/p|_initial
          NPT=0
          D0=YZXB(11)
        ENDIF

C----- File type zgoubi.spn
      ELSEIF(NL .EQ. NSPN) THEN   

        IF(IPASS.EQ. 1) NPT=0

      ENDIF  

      X = YZXB(KX) 
      Y = YZXB(KY) 

      X = X -DX * YZXB(1)
      Y = Y -DY * YZXB(1)

      IF(OKXAV) X=X-XMOY
      IF(OKYAV) Y=Y-YMOY
      IF(OKABSX) X = ABS(X)
      IF(OKABSY) Y = ABS(Y)
      IF(OKDX) THEN 
        BB = BX
        IF(READBX) BB = BTAB(ITAB,1)
        X = AX*(X-BB)**PX
      ELSEIF(OKABSX) THEN 
        BB = BX
        IF(READBX) BB = BTAB(ITAB,1)
        X = AX*X**PX + BB
      ELSEIF(OKABDX) THEN 
        BB = BX
        IF(READBX) BB = BTAB(ITAB,1)
        X = AX*ABS(X-BB)**PX
      ELSE
        BB = BX
        IF(READBX) BB = BTAB(ITAB,1)
        X = AX*X**PX + BB
      ENDIF
      IF(OKDY) THEN 
        BB = BY
        IF(READBY) BB = BTAB(ITAB,2)
        Y = AY*(Y-BB)**PY
      ELSEIF(OKABSY) THEN 
        BB = BY
        IF(READBY) BB = BTAB(ITAB,2)
        Y = AY*Y**PY + BB
      ELSEIF(OKABDY) THEN 
        BB = BY
        IF(READBY) BB = BTAB(ITAB,2)
        Y = AY*ABS(Y-BB)**PY
      ELSE
        BB = BY
        IF(READBY) BB = BTAB(ITAB,2)
        Y = AY*Y**PY + BB
      ENDIF
      
      ITAB = ITAB + 1

      IF    (NL .EQ. NPLT
     >           .OR. NL .EQ. NFAI) THEN
C------- Reset S at each pass
        IF(OKS) THEN
          IF(KX .EQ. 6) THEN
            X = X - SMAX*(IPASS-1)
C----------------
C Rustine LHC
C             IF(X .LT. 0.D0) X = 0.D0
C--------------------------
          ENDIF 
        ENDIF        
      ENDIF        

C----- File type zgoubi.plt
      IF    (NL .EQ. NPLT) THEN
          IF(NPT .NE. 0) THEN
            IF(X .LT. X0) NPT = 0
          ENDIF        
      ENDIF        

      IF( INRANG(X,Y,XMI,XMA,YMI,YMA) ) THEN

        NOC = NOC+1

        IF(NPT .GT. 0) THEN

          NPT = NPT+1
          YDX =  (X-X0) * Y
          SYDX = SYDX + YDX
          X0=X
            
          CALL VECTPL(X,Y,2) 

          IF(LIS .EQ. 2) THEN 
C            CALL IMPV(NLOG,NPT,X,Y,YDX,SYDX,NDX(2),NDX(1),KX,KY)
            CALL IMPV(NLOG,NOC,X,Y,YDX,SYDX,NDX(2))
          ENDIF

          IF(OKTAG) THEN
            IF(OKT(NDX(2)) .EQ. 0) THEN
              IF(X .GT. XMI+TAGX*(XMA-XMI)) THEN
                OKT(NDX(2)) = 1
                WRITE(TXT,FMT='(A)') LET 
                CALL TRTXT(X,Y+TAGY*(YMA - YMI),TXT,1)
              ENDIF
            ENDIF
          ENDIF

          XYLAST(1,NDX(2)) = X
          XYLAST(2,NDX(2)) = Y
          XM = XM + X
          YM = YM + Y
          X2M = X2M + X*X
          Y2M = Y2M + Y*Y
          XYM = XYM + X*Y

        ELSEIF(NPT .EQ. 0) THEN
C--------- Beginning of next optical lmnt or next pass
       
          NPT = 1

          YDX = 0.2D0
          SYDX = 0.D0
          X0 = X

            IF(CNECT) THEN
              IF(XYLAST(1,NDX(2)) .LE. X) THEN
C             ... protection ungainst previous traxe
                CALL VECTPL( XYLAST(1,NDX(2)), XYLAST(2,NDX(2)), 4)
              ELSE
                CALL VECTPL(X,Y,4)
              ENDIF
            ELSE
              CALL VECTPL(X,Y,4) 
            ENDIF

            CALL VECTPL(X,Y,2) 
            IF(LIS .EQ. 2)
     >            CALL IMPV(NLOG,NOC,X,Y,YDX,SYDX,NDX(2))
C     >            CALL IMPV(NLOG,NPT,X,Y,YDX,SYDX,NDX(2),NDX(1),KX,KY)       

            XYLAST(1,NDX(2)) = X
            XYLAST(2,NDX(2)) = Y
            XM = XM + X
            YM = YM + Y
            X2M = X2M + X*X
            Y2M = Y2M + Y*Y
            XYM = XYM + X*Y
        ENDIF

        IF(KLIPS) CALL FILCOO(KPS,NOC,YZXB)

      ENDIF

      IPASS0 = IPASS

      GOTO 44             
C     ------------------------------------

 19   CONTINUE
      CALL FBGTXT
      WRITE(6,*) ' SBR PLOTER: ERROR DURING READ OF EVENT #',NRD

 10   CONTINUE
      NPTR = NRD
      NPTS = NOC
      IF(NPTS.GT.NTRMAX) NPTS=NTRMAX
      CALL VECTPL(X,Y,4)
      IF(IT0 .NE. 0) THEN
        IF(KX.EQ. 7.AND.KY.EQ. 2..AND.PPA.NE.0.D0)
     >           CALL PPAR(PPA,X,Y,T0,KART)
      ENDIF

      CALL TRKVAR(NOC,KVAR(KY),KDIM(KY),KVAR(KX),KDIM(KX)) 
      WRITE(6,FMT='(/,'' SUM(Y)dX [XMI->XMA] ='',1P,E16.8,/)') SYDX 
      WRITE(6,*) ' PLOTTED ',NOC,' POINTS,  OVER',NRD-1,' READ IN FILE'

      AN = DBLE(NOC)
      XM = XM / AN
      YM = YM / AN
      X2M = X2M / AN
      Y2M = Y2M / AN
      XYM = XYM / AN

      SX = SQRT(X2M - XM*XM)
      SY = SQRT(Y2M - YM*YM)
      COV = XYM - XM*YM
      COR = COV / SX / SY

      IUNIT = 6
      WRITE(IUNIT,*) ' '
      WRITE(IUNIT,*) ' Means : <X>, SigmaX, <Y>, SigmaY, Correlation'
      WRITE(IUNIT,FMT='(1P,5G14.6)') XM,SX,YM,SY,COR
      WRITE(IUNIT,FMT=
     >'(''Linear regression :   Y ='',1P,E16.8,'' + '',E16.8,'' * X'')') 
     >YM-COR*SY/SX*XM, COR*SY/SX
      WRITE(IUNIT,*) ' '
C      IUNIT = NLOG
C      WRITE(IUNIT,FMT='(1P,5(1X,G14.6),1X,I7)') XM,SX,YM,SY,COR,IPASS
C      CALL FLUSH2(IUNIT,.FALSE.)      

      IF(LIS .EQ. 2) CLOSE(JUN)

      RETURN

      ENTRY PLOTE1
 11     WRITE(6,117) TAGX, TAGY
 117    FORMAT('   Give position coefficients for curves'' tags ',
     >  /,'         (normally,   0 < coeff_X < 1   &   |coeff_Y| < .1)',
     >  /,'   Now, in X & Y respectively : ',F10.4,'    &  ',F10.4,
     >  /,'   - Entering negative value for X  will cancel tagging -')
        READ(5,FMT='(2E12.4)',ERR=11) TAGX,TAGY
        IF(TAGX .GT. 1. .OR. ABS(TAGY) .GT. .1) GOTO 11
        OKTAG = TAGX .GE. 0.
      RETURN

      ENTRY PLOTE2
 20     WRITE(6,120) PPA
 120    FORMAT('  Give length for path extrapolation (cm)',
     >    '- Now :',1P,G12.4)
        READ(5,*,ERR=20) PPA
      RETURN

      ENTRY PLOTE3
        OKS0 = OKS
        WRITE(6,130)
 130    FORMAT('  Reset S coordinate to 0 at each pass (n/y)')
        READ(5,*,ERR=20) TXT
        OKS = TXT.EQ.'Y' .OR. TXT.EQ.'y'
        IF(OKS .NEQV. OKS0) OKECH = .FALSE.
      RETURN

      ENTRY PLOT31
        OKS=.FALSE.
      RETURN

      ENTRY PLOT4(OKXAVI)
        OKXAV=OKXAVI
      RETURN
      ENTRY PLOT41(OKYAVI)
        OKYAV=OKYAVI
      RETURN
      ENTRY PLOT4R(
     >             OKXAVO,OKYAVO)
        OKXAVO=OKXAV
        OKYAVO=OKYAV
      RETURN

      ENTRY PLOT5(AXI,PXI,BXI,DXI,AYI,PYI,BYI,DYI,IOPT,IBXY,IUN)
        AX = AXI
        BX = BXI
        PX = PXI
        DX = DXI
        AY = AYI
        BY = BYI
        PY = PYI
        DY = DYI
        OKDX = IOPT.EQ.2
        OKABSX = IOPT.EQ.3
        OKABDX = IOPT.EQ.4
        OKDY = IOPT.EQ.6
        OKABSY = IOPT.EQ.7
        OKABDY = IOPT.EQ.8
C        OKDX = OKDX .OR. IOPT.EQ.2
C        OKABSX = OKABSX .OR. IOPT.EQ.3
C        OKABDX = OKABDX .OR. IOPT.EQ.4
C        OKDY = OKDY .OR. IOPT.EQ.6
C        OKABSY = OKABSY .OR. IOPT.EQ.7
C        OKABDY = OKABDY .OR. IOPT.EQ.8
        IF(IBXY.NE.0) THEN
          READBX = IBXY.EQ.1
          READBY = IBXY.EQ.2
          ISTEP=1

 51       CONTINUE
C Units in BTAB should be m, as delivered by READCO
            READ(IUN,*,ERR=519,END=518) DUM,BTAB(ISTEP,IBXY)
            ISTEP=ISTEP+1
            GOTO 51
 518        WRITE(6,*) ' SBR PLOTER : end read BTAB upon EOF'
            GOTO 510
 519        WRITE(6,*) ' SBR PLOTER : end read BTAB upon error'
 510        CONTINUE
            WRITE(6,*) '  ->  ',ISTEP-1,' data have been read'
        ENDIF
      RETURN

      ENTRY PLOT51(AXJ,PXJ,BXJ,DXJ,AYJ,PYJ,BYJ,DYJ)
        AX = AXJ
        BX = BXJ
        PX = PXJ
        DX = DXJ
        AY = AYJ
        BY = BYJ
        PY = PYJ
        DY = DYJ
        OKDX = .FALSE.
        OKABSX = .FALSE.
          READBX = .FALSE.
        OKDY = .FALSE.
        OKABSY = .FALSE.
          READBY = .FALSE.
      RETURN

      ENTRY PLOT6(KLIPSI)
      KLIPS=KLIPSI
      RETURN

      ENTRY PLOT7(OKX12I)
        OKX12=OKX12I
      RETURN
      ENTRY PLOT71(OKY12I)
        OKY12=OKY12I
      RETURN
      ENTRY PLOT7R(
     >             OKX12O,OKY12O)
        OKX12O=OKX12
        OKY12O=OKY12
      RETURN

      ENTRY PLOT9(
     >            JPASS)
        JPASS = IPASS
      RETURN

C----- Scale computer
      ENTRY CALECH(NL,
     >                NOCE,NPASS)
      IF(KY.EQ. 28) THEN
C------- Histogram
        IF(.NOT.OKBIN) CALL BIN(NL,OKECH,KX,NB,
     >                           NOCE,OKBIN,XMOY,SIG,XMI,XMA,*98)
        YMI = 1.D10
        YMA =-1.D10
        DO 111 I=1,NB
          Y=BINS(I)
          IF(Y .GT. YMA) YMA=Y
          IF(Y .LT. YMI) YMI=Y
 111    CONTINUE
        CALL BIN3R(BINM)
        IF(.NOT.BINM) YMI=0.D0
        XMI=BINS(MXB+1)
        XMA=XMI + (NB-1.D0)*BINS(MXB+2)

        IF(OKXAV) THEN
          DO I=1,NB
            BINS(NB) = BINS(NB) - XMOY
          ENDDO
          BINS(MXB+1) = BINS(MXB+1) - XMOY
        ENDIF

      ELSE
        NOCE = 0
        XMOY=0.D0
        YMOY=0.D0
        XMI = 1.D10
        XMA =-1.D10
        YMI = 1.D10
        YMA =-1.D10
        NRD=0
        ITAB = 1
        CALL REWIN2(NL,*97)

C       ******** BOUCLE SUR READ FICHIER NL 
 144    CONTINUE                 
        NRD=NRD+1
        CALL READCO(NL,
     >                 KART,LET,YZXB,NDX,*110,*96) 
C        IF(NDX(1).LT. -1) GOTO 144

        IPASS=NINT(YZXB(39))
        NOCE = NOCE+1

      IF(IPASS .NE. IPASS0) NPT = 0

C----- File type zgoubi.plt
      IF    (NL .EQ. NPLT) THEN

C-------  Particle changes
        IF(LET .NE. LET0 .OR. NDX(2) .NE. IT0 ) THEN

C--------- If this is not the first particle...
C          IF(LET0 .NE. ' ') THEN
          IF(IT0 .NE. 0) THEN

C----------- Extrapolate the plot linearly over an additive path length PPA
C            at end of field maps
            IF(KX.EQ. 8 .AND. KY.EQ. 2..AND.PPA.NE. 0.D0)
     >        CALL PPAR(PPA,X,Y,T0,KART)
          ENDIF

          NPT=0
          itab = 1
          LET0=LET
          IT0 = NDX(2)
           
        ENDIF
        T0 = YZXB(3)

C----- File type zgoubi.fai
      ELSEIF(NL .EQ. NFAI) THEN   

        IF(YZXB(11) .NE. D0) THEN
C---------- Detects a change of dp/p|_initial
          NPT=0
          D0=YZXB(11)
        ENDIF

C----- File type zgoubi.spn
      ELSEIF(NL .EQ. NSPN) THEN   

        IF(IPASS.EQ. 1) NPT=0

      ENDIF  

        X = YZXB(KX)
        Y = YZXB(KY)

      X = X -DX * YZXB(1)
      Y = Y -DY * YZXB(1)

        XMOY = XMOY+X
        YMOY = YMOY+Y

      IF(OKABSX) X = ABS(X)
      IF(OKABSY) Y = ABS(Y)
      IF(OKDX) THEN 
        BB = BX
        IF(READBX) BB = BTAB(ITAB,1)
        X = AX*(X-BB)**PX
      ELSEIF(OKABSX) THEN 
        BB = BX
        IF(READBX) BB = BTAB(ITAB,1)
        X = AX*X**PX + BB
      ELSEIF(OKABDX) THEN 
        BB = BX
        IF(READBX) BB = BTAB(ITAB,1)
        X = AX*ABS(X-BB)**PX
      ELSE
        BB = BX
        IF(READBX) BB = BTAB(ITAB,1)
        X = AX*X**PX + BB
      ENDIF
      IF(OKDY) THEN 
        BB = BY
        IF(READBY) BB = BTAB(ITAB,2)
        Y = AY*(Y-BB)**PY
      ELSEIF(OKABSY) THEN 
        BB = BY
        IF(READBY) BB = BTAB(ITAB,2)
        Y = AY*Y**PY + BB
      ELSEIF(OKABDY) THEN 
        BB = BY
        IF(READBY) BB = BTAB(ITAB,2)
        Y = AY*ABS(Y-BB)**PY
      ELSE
        BB = BY
        IF(READBY) BB = BTAB(ITAB,2)
        Y = AY*Y**PY + BB
      ENDIF
      ITAB = ITAB + 1

c        X = AX*X**PX + BX
c        Y = AY*Y**PY + BY
        CALL MINMAX(X,Y,XMI,XMA,YMI,YMA)

        GOTO 144             

C       -----------------------------------
 110    CONTINUE

        NPASS = NINT(YZXB(39))

        XMOY = XMOY/NOCE
        YMOY = YMOY/NOCE
         
      ENDIF

      IF(OKXAV) THEN
        XMI=XMI-XMOY
        XMA=XMA-XMOY
      ENDIF
      IF(OKYAV) THEN
        IF(KY.NE. 28) THEN
          YMI=YMI-YMOY
          YMA=YMA-YMOY
        ENDIF
      ENDIF

C----- Sets SMAX for possible reset of S at each Pass (Menu-7/3/11)
      IF(KX .EQ. 6) THEN
C------- YZXB(39) is current value of IPASS 
        SMAX = XMA/YZXB(39)   
CCC        SMAX = (XMA-XMI)/YZXB(20)   
        DX = (SMAX - XMI) / 100.D0
        IF(OKS) THEN
          XMA = SMAX + DX
CCC          XMA = XMI + SMAX + DX
          XMI = XMI - DX
        ENDIF
      ENDIF
C      CALL TRKVA3(XMI,XMA,YMI,YMA)
      RETURN
                
 98   WRITE(6,*) ' *** Error Xmi-max while binning histo'
      RETURN
 96   WRITE(6,*) ' SBR CALECH, at event #',NRD,' :'
 97   WRITE(6,*) '     *** Error with data file (check existence ?)'
 99   RETURN 
      END
