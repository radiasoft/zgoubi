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
      SUBROUTINE BIN(NL,OKECH,KX,NB,
     >                               NOC,OKBIN,XMOY,SIG,XMI,XMA,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKBIN, OKECH
C     --------------------------  
C     REMPLI LES BINS POUR HISTO
C     --------------------------  
      PARAMETER (MXB=1000)
      COMMON/B/BINS(MXB+2)
      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*7, KPOL(2)*9, KDIM(MXVAR)*7
      COMMON/INPVR/ KVAR, KPOL, KDIM
      COMMON/REBELO/ NRBLT,IPASS,KWRI,NNDES,STDVM

      CHARACTER LET
      DIMENSION YZXB(MXVAR),NDX(5)
      SAVE NRD,XM,XM2

      LOGICAL BINM,BINMO
      SAVE BINM

      DIMENSION NOCY(MXB)

      INCLUDE "MXSTEP.H"
      DIMENSION BTAB(MXSTEP,2)
      SAVE AX,PX,BX,AY,PY,BY, BTAB
      DATA AX,PX,BX,AY,PY,BY / 1.D0, 1.D0, 0.D0, 1.D0, 1.D0, 0.D0 /

      LOGICAL  OKDX, OKABSX, OKABDX, OKDY, OKABSY, OKABDY
      SAVE     OKDX, OKABSX, OKABDX, OKDY, OKABSY, OKABDY
      DATA     OKDX, OKABSX, OKABDX, OKDY, OKABSY, OKABDY / 6*.FALSE. /

      LOGICAL READBX, READBY
      SAVE READBX, READBY
      DATA READBX, READBY / 2*.FALSE. /

      OKBIN=.FALSE.
      WRITE(6,*) ' Bining  for  histogram, wait ...'

      IF(.NOT. OKECH) THEN
         NRD=0
C------- Searching XMI, XMAX
          XMI = 1.D10
          XMA =-1.D10

          CALL REWIN2(NL,*96)
          NOC=0
          NRD=0
 
C--------- Loop over NL file
 45       CONTINUE
          NRD=NRD+1                 
          CALL READCO(NL,
     >                      KART,LET,YZXB,NDX,*11,*89) 
          IPASS=YZXB(39)
          NOC=NOC+1

          XX = YZXB(KX)
          IF(OKABSX) XX = ABS(XX)
          IF(OKDX) THEN 
            BB = BX
            IF(READBX) BB = BTAB(ITAB,1)
            XX = AX*(XX-BB)**PX
          ELSEIF(OKABSX) THEN 
            BB = BX
            IF(READBX) BB = BTAB(ITAB,1)
            XX = AX*XX**PX + BB
          ELSEIF(OKABDX) THEN 
            BB = BX
            IF(READBX) BB = BTAB(ITAB,1)
            XX = AX*ABS(XX-BB)**PX
          ELSE
            BB = BX
            IF(READBX) BB = BTAB(ITAB,1)
            XX = AX*XX**PX + BB
          ENDIF
          CALL MINMAX(XX,0.D0,XMI,XMA,DUM,DUM)
          GOTO 45             
C-----------------------------------------------

 89       WRITE(6,*) ' *** SBR BIN: Error during read of event #',NRD
 11       CONTINUE

      ENDIF  ! OKECH

      CALL RAZ(BINS,NB)
      IF (XMA.LE. XMI) GOTO 97
      BINS(MXB+1)=XMI
      DB=(XMA-XMI)/FLOAT(NB-1)
      BINS(MXB+2)=DB

C----- Fills up  the histogram
      CALL REWIN2(NL,*96)
      NOC=0
      NRD=0
      XM=0.D0
      XM2=0.D0
      IBMI = MXB+1
      IBMA = -1
      CALL IRAZ(NOCY,NB)

C----- Loop  over  READ file NL
 44   CONTINUE                 
      NRD=NRD+1
      CALL READCO(NL,
     >                  KART,LET,YZXB,NDX,*10,*10)
C      IPASS=NINT(YZXB(20))
      IPASS=NINT(YZXB(39))
      NOC=NOC+1
          XX = YZXB(KX)
          IF(OKABSX) XX = ABS(XX)
          IF(OKDX) THEN 
            BB = BX
            IF(READBX) BB = BTAB(ITAB,1)
            XX = AX*(XX-BB)**PX
          ELSEIF(OKABSX) THEN 
            BB = BX
            IF(READBX) BB = BTAB(ITAB,1)
            XX = AX*XX**PX + BB
          ELSEIF(OKABDX) THEN 
            BB = BX
            IF(READBX) BB = BTAB(ITAB,1)
            XX = AX*ABS(XX-BB)**PX
          ELSE
            BB = BX
            IF(READBX) BB = BTAB(ITAB,1)
            XX = AX*XX**PX + BB
          ENDIF
      IB=1+INT((XX-XMI)/DB)
C      IB=1+INT((YZXB(KX)-XMI)/DB)
      IF(IB.GT.IBMA) IBMA=IB
      IF(IB.LT.IBMI) IBMI=IB
      IF(BINM) THEN
C------- Bin average value of KY
          YY = YZXB(KYB)
          IF(OKABSY) YY = ABS(YY)
          IF(OKDY) THEN 
            BB = BY
            IF(READBY) BB = BTAB(ITAB,1)
            YY = AY*(YY-BB)**PY
          ELSEIF(OKABSY) THEN 
            BB = BY
            IF(READBY) BB = BTAB(ITAB,1)
            YY = AY*YY**PY + BB
          ELSEIF(OKABDY) THEN 
            BB = BY
            IF(READBY) BB = BTAB(ITAB,1)
            YY = AY*ABS(YY-BB)**PY
          ELSE
            BB = BY
            IF(READBY) BB = BTAB(ITAB,1)
            YY = AY*YY**PY + BB
          ENDIF
        DY =  YY
        NOCY(IB) = NOCY(IB)+1.D0
      ELSE
C------- Just COUNT events
        DY = 1.D0
        DY = AY * DY**PY + BY
      ENDIF
      BINS(IB)=BINS(IB) + DY
      XM=XM+XX
C      XM=XM+YZXB(KX)
      XM2=XM2+XX*XX
C      XM2=XM2+YZXB(KX)*YZXB(KX)
      GOTO 44             
C-----------------------------------------

 10   CONTINUE

      OKBIN=.TRUE.
      IF(BINM) THEN
        DO 9 IB = IBMI,IBMA
          IF(NOCY(IB).NE.0) BINS(IB)=BINS(IB)/NOCY(IB)
 9      CONTINUE
      ENDIF

      ENTRY BIN2(NOC,
     >               XMOY,SIG)

      XMOY = XM / NOC
      SIG=SQRT(XM2/NOC-XMOY**2)

      IF(NOC .EQ. NRD-1) THEN
        WRITE(6,*) NOC,' COUNTS HAVE BEEN BINED'
      ELSE
        WRITE(6,*) NOC,' COUNTS HAVE BEEN BINED OVER ',NRD-1,' FILED'
        WRITE(6,*) '(',NRD-1-NOC,' COUNTS HAVE KEX<0 )'
      ENDIF
      RETURN

      ENTRY BIN3W(MBIN,KYBI) 
      BINM = MBIN .EQ. 1
      KYB = KYBI
      RETURN         
         
      ENTRY BIN3R(BINMO)
      BINMO = BINM
      RETURN         
         
      ENTRY BIN4W(AXI,PXI,BXI,AYI,PYI,BYI,IOPT,IBXY,IUN)
        AX = AXI
        BX = BXI
        PX = PXI
        AY = AYI
        BY = BYI
        PY = PYI
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
c            write(99,*) DUM,BTAB(ISTEP,IBXY), ' sbr ploter, test plot5'
            ISTEP=ISTEP+1
            GOTO 51
 518        WRITE(6,*) ' SBR PLOTER : end read BTAB upon EOF'
            GOTO 510
 519        WRITE(6,*) ' SBR PLOTER : end read BTAB upon error'
 510        CONTINUE
            WRITE(6,*) '  ->  ',ISTEP-1,' data have been read'
        ENDIF
      RETURN

      ENTRY BIN41W(AXJ,PXJ,BXJ,AYJ,PYJ,BYJ)
        AX = AXJ
        BX = BXJ
        PX = PXJ
        AY = AYJ
        BY = BYJ
        PY = PYJ
        OKDX = .FALSE.
        OKABSX = .FALSE.
          READBX = .FALSE.
        OKDY = .FALSE.
        OKABSY = .FALSE.
          READBY = .FALSE.
      RETURN

 97   WRITE(6,*) ' SBR BIN: ERROR   XMI > XMA'      
      WRITE(6,FMT='(/,I7,''  particles  counted'',/)') NOC
      RETURN 1
 96   WRITE(6,*) ' SBR BIN: Error during read in data file'
      END
