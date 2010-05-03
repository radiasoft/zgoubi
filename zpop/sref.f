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
      SUBROUTINE SREF(KEP,OX,IX,IY,Q,AM,FE,NRMA,NRMC,
     >                                           GAM,R,NOC,NRD,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------------------------------------------
C     Calculates E field from synchrotron radiation
C     Upon option (KEP) calculate scales, or plot
C     ---------------------------------------------
      DIMENSION OX(*)
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN
      COMMON/LUN/NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      INCLUDE 'MXSTEP.H'
      PARAMETER (MSAM=1000)
      COMMON/SRTAB/ EFD(MXSTEP,3),D3W(MSAM,9),NPTS
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      CHARACTER LET, TXT*80, REP
      INCLUDE 'MXVAR.H'
      DIMENSION YZXB(MXVAR),NDX(5)

      INCLUDE 'MAXNTR.H'
C      DOUBLE PRECISION SYDX

C-----          Part.   Part.-Obs.                           Normalized
C               posit.   direct      beta   beta-dot  field  obsrvr posit.
      DIMENSION RP(3),   RO(3),      BE(3),  BED(3),   B(3),  ON(3)

C----- Radiated E-field
C      in observer time
      DIMENSION E(3),  ONB(3), CX(3), CY(10)
      LOGICAL NEWL, NEWT

      SAVE MXSTP


      NOC = 0
      DCX = 0.D0

      IF(KEP .EQ. 1) CALL SRPLT(NLOG,KEP,DCX,CX,CY,IX,IY,UNB,RP,TPART,
     >                                               NOC,NPTS,NRMA,NRMC,
     >                                                   TM,YM,TC,SYDX)

C      S0 = -1.D10
      S0= 0.D0
      NOV = 0
      NRD = 0
      NOEL = 0
      NTRAJ = 0
      NEWL = .FALSE.

C----- BOUCLE SUR READ FICHIER NL 
 44   CONTINUE                 

        NRD=NRD+1
        IF(NRD .GT. MXSTP) GOTO 10

        CALL XYZBW(NRD,
     >                 KART,LET,YZXB,NDX)

        IF(NDX(1) .LT. -1) GOTO 95
        IF( KART .EQ. 2 ) GOTO 96
        NOC = NOC+1

 45     CONTINUE

        DS = YZXB(9)                     ! Size of next step (m)
        S = YZXB(6)                      ! Path length (m)
        B(1) = YZXB(30)                ! Bx, Tesla
        B(2) = YZXB(31)                ! By
        B(3) = YZXB(32)                ! Bz

        NEWL = NOEL .NE. NDX(5)
        NEWT = NTRAJ .NE. NDX(2) 

        IF(NEWT) THEN
          NTRAJ =  NDX(2)
          CALL SRPLTI
C           call fbgtxt
C              write(*,*)  ' -------  newt =',NEWT
        ENDIF

        IF(NEWL) THEN
          
          IF(NOC .GE. 2) THEN
C----------- Entering next magnetic element, account for a drift
C            between the two successive elements

C----------- Corrections
            EFD(NOC-1,1) = 0.D0 
            DCI = -CX(IX)

C----------- Begining of drift
C----------- Radiated electric field is zero
            E(1) = 0.D0
            E(2) = 0.D0
            E(3) = 0.D0

C----------- dt/dt' = 1 - n.beta
C----------- Account for a drift between the two successive elements
            DTI = (S - S0) / ( BTA * CL )                   ! Prtcl time step (s) 
            CX(1) = DTI * UNB                         ! dt at Observer (s) 

            CY(2) = E(2)                              ! Ey 
            CY(3) = E(3)                              ! Ez 

            IF( KEP .LE. 2 ) THEN
              CX(2) = DTI                             ! dt' at particle
              CY(1) = E(1)                           ! Ex 
              CALL SRPLT(NLOG,KEP,DCI,CX,CY,IX,IY,UNB,RP,TPART,
     >                                              NOC,NPTS,NRMA,NRMC,
     >                                                 TM,YM,TC,SYDX)
            ENDIF

C----------- Store time step dt at Observer
C            and  Ey, Ez for further Fourier transform
            IF ( NOC .LE. MXSTEP ) THEN
              EFD(NOC,1) = CX(1)                        ! dt at observer (s)
              EFD(NOC,2) = CY(2)                        ! Ey  (V/m)
              EFD(NOC,3) = CY(3)                        ! Ez  (V/m)
            ENDIF

C----------- End of drift

            NOC = NOC + 1
            DS = 0.D0
            B(1) = 0.D0
            B(2) = 0.D0
            B(3) = 0.D0

          ENDIF      
C--------------endif NOC

        ENDIF
C----------- endif NEWL

CC------- Test overlapping of fringe fields
C        IF((S - S0) .LT. 0.D0) CALL SRCUL(S0,
C     >                                 S,NOV,*98)
        S0 = S

        NOEL = NDX(5)

C------- Displacement step in magnet frame  ( dr(t') )
        CT = COS( YZXB(3) )
        CP = COS( YZXB(5) )
C        DX = DS * CT * CP

C------- Particle position in magnet frame  ( r(t') )
C        RP(1) = RP(1) + DX                        ! X (m)
C        RP(1) = YZXB(7)                   ! X (m)
        RP(1) = YZXB(8)                   ! X (m)
        RP(2) = YZXB(2)                   ! Y (m)
        RP(3) = YZXB(4)                   ! Z (m)
C             write(*,*)rp(1),rp(2),rp(3)
C------- Particle to observer vector  ( R(t') )
        RO(1) = OX(1) - RP(1)
        RO(2) = OX(2) - RP(2)
        RO(3) = OX(3) - RP(3)
        R = SQRT(RO(1)*RO(1) + RO(2)*RO(2) + RO(3)*RO(3))

C------- Normalized ( n(t') )
        ON(1) = RO(1)/R
        ON(2) = RO(2)/R
        ON(3) = RO(3)/R

C------- Beta
C   yzxb(1)=dp/p 
        BRO = YZXB(36) * (1.D0 + YZXB(1))        ! Rigidity, T.m
        P0 = BRO * CL * 1.D-6 * Q / QE           ! Momentum, MeV/c              
        BTA = P0 / SQRT( P0*P0 + AM*AM )  
C FM, LC 9/99
C        QMG = BTA * CL * CL / (P0 * 1.D6)        ! = q/m.gamma, MKSA units
        QMG = BTA * CL * CL / (P0 * 1.D6)* Q / QE ! = q/m.gamma, MKSA units

        BE(1) = BTA * CT * CP
        BE(2) = BTA * SIN( YZXB(3) ) * CP
        BE(3) = BTA * SIN( YZXB(5) )
                 
C------- Acceleration = Beta-dot = q/m beta x Field
        BED(1) = QMG * ( BE(2)*B(3) - BE(3)*B(2) )
        BED(2) = QMG * ( BE(3)*B(1) - BE(1)*B(3) )
        BED(3) = QMG * ( BE(1)*B(2) - BE(2)*B(1) )

C        BEn = SQRT(BE(1)*BE(1)+BE(2)*BE(2)+BE(3)*BE(3))
C        Bn = SQRT(B(1)*B(1)+B(2)*B(2)+B(3)*B(3))
C        BTAD = SQRT(BED(1)*BED(1)+BED(2)*BED(2)+BED(3)*BED(3))

C------- 1 - n.Beta
        UNB = (1.D0 - ON(1)*BE(1)) - (ON(2)*BE(2)+ON(3)*BE(3) )

        G1 = BTA * AM / P0
        GAM = 1.D0 / G1
CCCCCCCCCCC   Rustine superaco - interdit le passage 
C             Aussi, etude pb LHC a 30 TeV -> pics de bords dyssimetriques
C        IF(UNB*UNB .LT. 0.D0) THEN
CCCCCCCCCCC       
        IF(UNB*UNB .LT. 1.D-3) THEN

          PHI = ATAN( OX(2) / OX(1) )
          PSI = ATAN( OX(3) / SQRT( OX(1)*OX(1) + OX(2)*OX(2) ) )
          SPH = SIN( .5D0 * PHI )
          SPS = SIN( .5D0 * PSI )
          SPH2 = SPH * SPH
          SPS2 = SPS * SPS
          EPX = 2.D0 * (SPH2 + SPS2) - 4.D0 * SPH2 * SPS2
          G2 = G1 * G1 
          A = -( G2 + BE(2) * BE(2) + BE(3) * BE(3) )
          AN = .5D0
          FAC = 1.D0
          IQ = 1
          QSX = AN * A
 2        CONTINUE
            AN = AN * (AN - IQ)
            IQ = IQ + 1
            FAC = FAC * IQ
            DQ = AN / FAC * A**IQ  
            QSX = QSX + DQ
            TST = DQ/QSX
            TST = TST * TST
            IF(TST .LT. 1.D-16) GOTO 3
C            IF(TST .LT. 1.D-20) GOTO 3
            IF(IQ .GT. 100) GOTO 93
            GOTO 2
 3        CONTINUE
          QSX = -QSX
          UNB = ( EPX + QSX - ON(2)*BE(2) - ON(3)*BE(3) ) - EPX*QSX
          IF(UNB .LE. 0.D0) THEN
            CALL FBGTXT
            PRINT*,'  1 - n.beta  going negative at step NOC=',NOC
            print*,' unb, beta, R, on**2 ',
     >        unb, 
     >        sqrt(be(1)*be(1)+be(2)*be(2)+be(3)*be(3)),
     >             R,
     >        sqrt(on(1)*on(1)+on(2)*on(2)+on(3)*on(3))
C            print*,' epx, qsx, on_yz, be_yz ',
C     >        epx,qsx,on(2),be(2),on(3),be(3),unb
            PRINT*,' Try to increase observation distance (X,Y,Z)'
            GOTO 94
          ENDIF
          UNB2 = UNB * UNB
          UNB3 = UNB2 * UNB
C--------- n - Beta:
          ONB(1) = -EPX + QSX
          ONB(2) = ON(2) - BE(2)
          ONB(3) = ON(3) - BE(3)
C--------- (n - Beta) x Beta-dot
          V1 = ONB(2) * BED(3) - ONB(3) * BED(2)
          V2 = ONB(3) * BED(1) - ONB(1) * BED(3)
          V3 = ONB(1) * BED(2) - ONB(2) * BED(1)
        ELSE
          UNB2 = UNB * UNB
          UNB3 = UNB2 * UNB
C--------- (n - Beta) x Beta-dot
          V1 = (ON(2) - BE(2)) * BED(3) - (ON(3) - BE(3)) * BED(2)
          V2 = (ON(3) - BE(3)) * BED(1) - (ON(1) - BE(1)) * BED(3)
          V3 = (ON(1) - BE(1)) * BED(2) - (ON(2) - BE(2)) * BED(1)

        ENDIF

C------- Radiated electric field x,y,z components ( V/m )
        E(1) = FE * ( ON(2)*V3 - ON(3)*V2 ) / UNB3 / R
        E(2) = FE * ( ON(3)*V1 - ON(1)*V3 ) / UNB3 / R
        E(3) = FE * ( ON(1)*V2 - ON(2)*V1 ) / UNB3 / R
C              write(*,*) '   unb approx    ', UNB, 
C     >             sqrt(E(1)*E(1)+E(2)*E(2)+E(3)*E(3))
C              write(78,*) ' iq,  an, qsx ', iq, an,qsx, 
C            write(*,*) EPX,QSX,ON(2),BE(2),ON(3),BE(3),unb, 
C              write(*,*) ox(2),rp(1),rp(2),
C              write(*,*) on(2),ro(2),r,unb, 
C              write(*,*) ox(2),rp(2),ro(2),r,on(2), iq, 
C     >             sqrt(E(1)*E(1)+E(2)*E(2)+E(3)*E(3))

C------- dt/dt' = 1 - n.beta
        DTI = DS / ( BTA * CL )                   ! Prtcl time step dt' (s) 
        CX(1) = DTI * UNB                         ! dt at Observer (s) 
C        if(cx(1).le.0.d0) print*,ds,beta,dti,unb,cx(1)
        CY(2) = E(2)                              ! Ey 
        CY(3) = E(3)                              ! Ez 

        IF( KEP .LE. 2 ) THEN
C          CX(2) = S / ( BTA * CL )               ! Time t' at particle (s)
          CX(2) = DTI                            ! Prtcl time step dt' (s)
          CY(1) = E(1)                           ! Ex 

          CALL SRPLT(NLOG,KEP,DCX,CX,CY,IX,IY,UNB,RP,TPART,
     >                                    NOC,NPTS,NRMA,NRMC,
     >                                             TM,YM,TC,SYDX)
        ENDIF

C------- Store time step dt at Observer
C        and  Ey, Ez for further Fourier transform
        IF ( NOC .LE. MXSTEP ) THEN
          EFD(NOC,1) = CX(1)                        ! dt at observer (s)
          EFD(NOC,2) = CY(2)                        ! Ey  (V/m)
          EFD(NOC,3) = CY(3)                        ! Ez  (V/m)
        ENDIF

        IF(NEWL) GOTO 45

      GOTO 44             
C     ------------------------------------------

 10   CONTINUE

      CALL FBGTXT

      IF(KEP .EQ. 1) THEN

        IF(NRMA .EQ. 1) THEN

          WRITE(6,*) 'Normalization of scales:'
          WRITE(6,*)
          WRITE(6,*) '          --- X axis ---'
          WRITE(6,*) ' Calculated critical time, tc =',TC,' s'
          WRITE(6,*) ' ( omga_c = ',1.D0 / TC,' rad/s )'
          WRITE(6,*) ' Calculated time origin t(Ey_max): ',TM,' s'
 201      WRITE(6,*) ' Want to modify (N/Y) ?'
          READ(5,FMT='(A1)',ERR=201) REP 
          IF(REP .EQ. 'Y' .OR. REP .EQ. 'y') THEN
 202        WRITE(6,*) ' Give the value for time origin :'
            READ(5,*,ERR=202) TM
          ENDIF
 20       WRITE(6,*) ' Want to modify normalization tc (N/Y) ?'
          READ(5,FMT='(A1)',ERR=20) REP 
          IF(REP .EQ. 'Y' .OR. REP .EQ. 'y') THEN
 203        WRITE(6,*) ' Give the value for time normalisation tc:'
            READ(5,*,ERR=203) TC
          ENDIF

          IF(TC .NE. 0.D0) THEN
            XMI = (XMI-TM) / TC
            XMA = (XMA-TM) / TC
          ELSE
            IX = 1
            WRITE(6,*) ' Sorry, could not normalize the X-axis'
            WRITE(6,*) ' Will plot v.s. non-normalized observer time'
          ENDIF

          WRITE(6,*) '          --- Y axis ---'
          WRITE(6,*) ' Maximum field Ey_max = ',YM,' V/m'
 21       WRITE(6,*) ' Want to modify normalization Ey_max (N/Y) ?'
          READ(5,FMT='(A1)',ERR=21) REP 
          IF(REP .EQ. 'Y' .OR. REP .EQ. 'y') THEN
 19         WRITE(6,*) ' Give the value for normalization Ey_max :'
            READ(5,*,ERR=19) YM
          ENDIF
          IF(YM .NE. 0.D0) THEN
            YMI = YMI / YM
            YMA = YMA / YM
          ELSE
            IY = IY - 3
            WRITE(6,*) ' Sorry, could not normalize the Y-axis'
            WRITE(6,*) ' Will plot non-normalized E-field'
          ENDIF
         
          IF(NRMC .EQ. 1) NRMC = 0

        ELSEIF(NRMA .EQ. 2) THEN
        ENDIF

      ELSEIF(KEP .EQ. 2) THEN

        WRITE(TXT,FMT='('' SUM(Y)dX [XMI->XMA] ='',1P,E16.8)') SYDX
        CALL TRTXT(10.D0,10.D0,TXT,0)
        CALL FBGTXT
        WRITE(6,FMT='(/,'' SUM(Y)dX [XMI->XMA] ='',1P,E16.8,/)') SYDX
        WRITE(6,*) ' PLOTTED ',NOC,' POINTS, OVER ',
     >  NRD-1,' READ IN FILE'
        IF(NOV .GT. 0) 
     >  WRITE(6,*) ' ',NOV,' points overlap in fringe fields'
      ENDIF

      IF ( NOC .LE. MXSTEP ) THEN
        NPTS = NOC
        IF(NPTS.GT.NTRMAX) NPTS=NTRMAX
      ELSE 
        NPTS = MXSTEP
        IF(NPTS.GT.NTRMAX) NPTS=NTRMAX
        GOTO 97
      ENDIF
     
      RETURN
                 
 93   WRITE(6,*)
     >' SBR SREF: Process stopped; MORE THAN 100 iterations on QsiX'
      RETURN 1
 94   WRITE(6,*)
     >' SBR SREF: Process stopped; calculations go wrong (1-n.beta<0)'
      RETURN 1
 95   WRITE(6,*)
     >' SBR SREF: Process stopped; particle out of range (IEX=-1)'
      RETURN 1
 96   WRITE(6,*)
     >' SBR SREF: Process stopped; cannot deal with polar frames'
      RETURN 1
 97   WRITE(6,*) ' SBR SREF: Storage stopped at point ', MXSTEP
      WRITE(6,*) '           due to limited size of storage array'
      RETURN 1
C 98   WRITE(6,*)
C     >' SBR SREF: Process stopped; donnot overlap fringe fields'
C      RETURN 1

      ENTRY SREFW(MXSTPI)
      MXSTP = MXSTPI
      RETURN

      END
