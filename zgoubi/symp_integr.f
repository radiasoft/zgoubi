C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory 
C  C-AD, Bldg 911
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE SympIntegr(NRES, newIntegQ, MQI)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     (1) Optical elements defined in cartesian frame
C     (2) The drift-kick-drift symplectic particle motion integrator
C         for quadrupole (MQchoice = 1, newIntegQ = 1)
C     (3) The matrix-kick-matrix symplectic particle motion integrator
C         for quadrupole (MQchoice = 1, newIntegQ = 2)
C     (4) The drift-kick-drift symplectic particle motion integrator
C         for a general magnetic multipole (MQchoice = 2, newIntegQ = 1)
C-----------------------------------------------------------------------
      INCLUDE "MAXTRA.H"   ! PARAMETER (MXJ=7)
      INCLUDE "MAXCOO.H"   ! PARAMETER (MXT=10000)
      INCLUDE "MXFS.H"     ! PARAMETRE (MXF=45, MXS=2000, MLF=9, MXP=11)
      INCLUDE "MXLD.H"     ! PARAMETER (MXL=15000 , MXD=1400)
      INCLUDE "MXSCL.H"    ! PARAMETER (MXSCL=10)
      INCLUDE "MXSTEP.H"   ! PARAMETER (MXSTEP = 1000000)
      INCLUDE "C.CONST.H"  ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QE,AMPROT,CM2M
      INCLUDE "C.DON.H"    ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE 'C.SCAL.H'   ! COMMON/SCAL/ SCL(MXF,MXS,MXSCL),TIM(MXF,MXS),NTIM(MXF),KSCL
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"  ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),
                           ! IMAX,IEX(MXT),IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.PTICUL.H" ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE 'C.INTEG.H'  ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      INCLUDE "C.OBJET.H"  ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT,KZOB
      INCLUDE "C.RIGID.H"  ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE "C.UNITS.H"  ! COMMON/UNITS/ UNIT(MXJ)

      DOUBLE PRECISION XL, XLE, RO, GAP, ROBR, XS, YS, ZR, ZS, YR
      DOUBLE PRECISION DLE, XLS, DLS, DL0, TCUM, SCUM, l0, w0, eb2
      DOUBLE PRECISION FINTE, FINTS, PREF, EREF, BETA0, GAMMA0
      DOUBLE PRECISION b0g0, step, radius, mass, DUMMY, XLT, EXPd
      DOUBLE PRECISION PAS_E, PAS_L, PAS_S, step_E, step_L, step_S
      DOUBLE PRECISION Xpos, dist, G0E, G0S, G2E, G2S, G4E, G4S
      DOUBLE PRECISION Bfield(10), RSA(10), eb(10), XLtrc(MXT), ppii
      DOUBLE PRECISION dle1, dle2, dle3, dle4, dls1, dls2, dls3, dls4
      DOUBLE PRECISION Pd0, Pd1, Pd2, Pd3, Pd4, P1G1, P1G2, P2G2, P3G3
C      DOUBLE PRECISION EXPd
      DOUBLE PRECISION G0(MXSTEP), G2(MXSTEP), G4(MXSTEP), inv_beta0
      DOUBLE PRECISION eb2G0(MXSTEP), eb2G2(MXSTEP), eb2G4(MXSTEP)
      DOUBLE PRECISION FFEE(10), CE(0:5), FFES(10), CS(0:5)
      INTEGER KUASEX, NRES, newIntegQ, IT, NCE, NCS, MQI, HOMP
      INTEGER NSTEP, NSTEP_E, NSTEP_L, NSTEP_S, NSTEP_EL
      INTEGER J, JJ, JJ1, JJ2, JJ3, JJ4

      PARAMETER (l0 = 100.D0)
      PARAMETER (w0 = 299792458.D0)

C-------- Get the parameters for magnetic quadrupole (MQI .EQ. 1) ------
C---------or magnetic multipole (MQI .EQ. 2) ---------------------------
      IF (MQI .EQ. 1) THEN
        XL = A(NOEL,10)
        RO = A(NOEL,11)
        Bfield(2) = A(NOEL,12)

        XLE = A(NOEL,20)                     ! Integration zone: entrance
        DLE = A(NOEL,21)                     ! Dipole fringe field lambda_E
        NCE = NINT(A(NOEL,30))               ! NCE = number of coefs. for entrance fringe field
        XLS = A(NOEL,40)                     ! Integration zone: exit
        DLS = A(NOEL,41)                     ! Dipole fringe field lambda_S
        NCS = NINT(A(NOEL,50))               ! NCS = number of coefs. for entrance fringe field
        DO J = 0, 5
          CE(J) = A(NOEL,31+J)               ! Entrance fringe field coefficients
          CS(J) = A(NOEL,51+J)               ! Exit fringe field coefficients
        ENDDO

        PAS = A(NOEL,60)
        KP  = NINT(A(NOEL,70))
        XCE = A(NOEL,71)
        YCE = A(NOEL,72)
        ALE = A(NOEL,73)

        IF(NRES.GT.0) THEN
          WRITE(NRES, 100) XL, RO, Bfield(2), Bfield(2)/RO
        ENDIF
        KUASEX = 2
        HOMP = 2

      ELSE IF(MQI .EQ. 2) THEN
        XL = A(NOEL,2)                       ! Length
        RO = A(NOEL,3)                       ! Radius

        DO J = 1, 10
          Bfield(J) = A(NOEL,3+J)            ! Magnetic field from 2J-pole
          RSA(J) = A(NOEL,49+J)              ! Skew angles of field components
          IF (J .GT. 1) THEN
            FFEE(J) = A(NOEL,14+J)           ! entrance fringe field extend = lambda_E * FFE(J != 1)
            FFES(J) = A(NOEL,32+J)           ! exit fringe field extend = lambda_S * FFS(J != 1)
          ENDIF
          IF (J .LE. 6) THEN
            CE(J-1) = A(NOEL,25+J)           ! Entrance fringe field coefficients
            CS(J-1) = A(NOEL,43+J)           ! Exit fringe field coefficients
          ENDIF
        ENDDO

        XLE = A(NOEL,14)                     ! Integration zone: entrance
        DLE = A(NOEL,15)                     ! Dipole fringe field lambda_E
        NCE = NINT(A(NOEL,25))               ! NCE = number of coefs. for entrance fringe field
        XLS = A(NOEL,32)                     ! Integration zone: exit
        DLS = A(NOEL,33)                     ! Dipole fringe field lambda_S
        NCS = NINT(A(NOEL,43))               ! NCS = number of coefs. for exit fringe field

        PAS = A(NOEL,60)                     ! Integration step
        KP = NINT(A(NOEL,63))                ! In the Zgoubi manu, KPOS is located at # 61
                                             ! but in the program, it was move by two spots
        IF ((KP .EQ. 1) .OR. (KP .EQ. 2) .OR. (KP .EQ. 3)) THEN
          XCE = A(NOEL,64)
          YCE = A(NOEL,65)
          ALE = A(NOEL,66)
        ELSE IF (KP. EQ. 4) THEN
          XS = A(NOEL, 64)
          YS = A(NOEL, 65)
          ZR = A(NOEL, 66)
          ZS = A(NOEL, 67)
          YR = A(NOEL, 68)
        ELSE
          CALL ENDJOB('The value of KPOS is invalid', -99)
        ENDIF

        IF(NRES.GT.0) THEN
          WRITE(NRES, 101) XL, RO, (Bfield(J), J=1, 10)
        ENDIF

        KUASEX = 0                           ! Lowest-order multipole
        HOMP = 0                             ! Highest-order multipole
        J = 1
        DO WHILE (J .LE. 10)
          IF (Bfield(J) .NE. 0.D0) THEN
            IF (KUASEX .EQ. 0) KUASEX = J
            HOMP = J
          ENDIF
          J = J + 1
        END DO
        WRITE (*,99) KUASEX, HOMP

      ELSE
        CALL ENDJOB('MQI cannot be this value', -99)
      ENDIF

      IF (KUASEX .NE. 0) THEN
        GAP = RO/KUASEX
      ELSE
        CALL ENDJOB(' KUASEX = 0, all field is zero', -99)
      ENDIF

C      NSTEP_E = 0
C      NSTEP_S = 0
      NSTEP_L = CEILING(XL/PAS)
      IF (NSTEP_L .GT. MXSTEP) THEN
        CALL ENDJOB('Number of integration step > MASTEP', -99)
      ENDIF
      PAS_L = XL/NSTEP_L
      NSTEP = NSTEP_L

      IF((DLE .LT. 0.D0) .OR. (DLS .LT. 0.D0)) CALL
     >  ENDJOB('Delta_E or Delta_S cannot be negative', -99)
      IF((XLE .LT. 0.D0) .OR. (XLS .LT. 0.D0)) CALL
     >  ENDJOB('XLE or XLS cannot be negative', -99)

      DL0 = DLE + DLS
      IF(DL0 .EQ. 0.D0) THEN
C-------- Sharp edge at entrance and exit-------------------------------
        FINTE = XLE
        XLE = 0.D0
        FINTS = XLS
        XLS = 0.D0
        IF(NRES.GT.0) THEN
          WRITE(NRES,120) FINTE, FINTS, GAP
        ENDIF
      ELSE
C-------- Fringe fields are present-------------------------------------
        IF (MQI .EQ. 1) THEN            ! Quadrupole
          IF(NRES.GT.0) THEN
            IF((DLE .EQ. 0.D0) .OR. (XLE .EQ. 0.D0)) THEN
              FINTE = XLE
              XLE = 0.D0
              WRITE(NRES,121) FINTE, GAP
            ELSE
              WRITE(NRES,128)
              WRITE(NRES,130) XLE, DLE, NCE, (CE(J), J=0, NCE-1)
              WRITE(NRES,131)
C              NSTEP_E = CEILING (XLE/PAS)
C              PAS_E = XLE/NSTEP_E       ! Integration step for entrance fringe field
C              step_E = PAS_E/l0         ! Scaled integration step for entrance fringe field
C              dle1 = DLE/l0
C              dle2 = dle1*dle1
C              dle3 = dle2*dle1
C              dle4 = dle3*dle1
            ENDIF

            IF((DLS .EQ. 0.D0) .OR. (XLS .EQ. 0.D0)) THEN
              FINTS = XLS
              XLS = 0.D0
              WRITE(NRES,122) FINTS, GAP
            ELSE
              WRITE(NRES,129)
              WRITE(NRES,130) XLS,DLS, NCE, (CS(J), J=0, NCS-1)
              WRITE(NRES,131)
C              NSTEP_S = CEILING (XLS/PAS)
C              PAS_S = XLS/NSTEP_S       ! Integration step for exit fringe field
C              step_S = PAS_S/l0         ! Scaled integration step for exit fringe field
C              dls1 = DLS/l0
C              dls2 = dls1*dls1
C              dls3 = dls2*dls1
C              dls4 = dls3*dls1
            ENDIF
          ENDIF

C          NSTEP_EL = NSTEP_E + NSTEP_L
C          NSTEP = NSTEP_EL + NSTEP_S
          IF (NSTEP .GT. MXSTEP) THEN
            CALL ENDJOB('Number of integration step > MASTEP', -99)
          ENDIF

C          DO J = 1, NSTEP
C            IF (J .LE. NSTEP_E) THEN              ! Entrance fringe field
C              Xpos = PAS_E*(J - NSTEP_E - 0.5D0)
C            ELSE IF (J .GT. NSTEP_EL) THEN        ! Exit fringe field
C              Xpos = XL + PAS_S*(J - NSTEP_EL - 0.5D0)
C            ELSE
C              Xpos = PAS_L * (J - NSTEP_E - 0.5D0)
C            ENDIF

C            IF ((ABS(Xpos) .LE. XLE) .AND. (DLE .GT. 0.D0)) THEN
C              dist = -Xpos/DLE
C              Pd0 = 0.D0
C              Pd1 = 0.D0
C              Pd2 = 0.D0
C              Pd3 = 0.D0
C              Pd4 = 0.D0
C              DO JJ = 0, NCE-1
C                JJ1 = JJ
C                JJ2 = JJ1*(JJ-1)
C                JJ3 = JJ2*(JJ-2)
C                JJ4 = JJ3*(JJ-3)
C                Pd0 = Pd0 + CE(JJ)*(dist**JJ)
C                Pd1 = Pd1 + JJ1*CE(JJ)*(dist**(JJ-1))/dle1
C                Pd2 = Pd2 + JJ2*CE(JJ)*(dist**(JJ-2))/dle2
C                Pd3 = Pd3 + JJ3*CE(JJ)*(dist**(JJ-3))/dle3
C                Pd4 = Pd4 + JJ4*CE(JJ)*(dist**(JJ-4))/dle4
C              ENDDO
C
C              IF (Pd0 .GT. 40.D0) THEN           ! Avoid IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
C                G0E = 0.D0                       ! for generating a number < 1.D-38
C                G2E = 0.D0                       ! Second-order derivative for the FF
C                G4E = 0.D0                       ! Fourth-order derivative
C              ELSE IF (Pd0 .LT. -40.D0) THEN
C                G0E = 1.D0
C                G2E = 0.D0
C                G4E = 0.D0
C              ELSE
C                EXPd = EXP(Pd0)
C                G0E  = 1.D0/(1.D0 + EXPd)
C                P1G1 = EXPd*G0E                  ! EXPd*G0E
C                P1G2 = P1G1*G0E                  ! EXPd*(G0E**2)
C                P2G2 = P1G2*EXPd                 ! (EXPd**2)*(G0E**2)
C                P3G3 = P1G1*P2G2                 ! (EXPd**3)*(G0E**3)
C
C                G2E = P1G2*((-1.D0+2.D0*P1G1)*(Pd1**2) - Pd2)
C                G4E = P1G2*((-1.D0+14.D0*P1G1-36.D0*P2G2+24.D0*P3G3)
C     > *(Pd1**4) + (-6.D0+36.D0*P1G1-36.D0*P2G2)*(Pd1**2)*Pd2
C     > + (-4.D0+8.D0*P1G1)*Pd1*Pd3 + (-3.D0+6.D0*P1G1)*(Pd2**2)
C     > - 2.D0*Pd4)
C              ENDIF
C            ELSE
C              G0E = 1.D0
C              G2E = 0.D0
C              G4E = 0.D0
C            ENDIF
C
C            IF ((ABS(Xpos-XL) .LE. XLS) .AND. (DLS .GT. 0.D0)) THEN
C              dist = (Xpos-XL)/DLS
C              Pd0 = 0.D0
C              Pd1 = 0.D0
C              Pd2 = 0.D0
C              Pd3 = 0.D0
C              Pd4 = 0.D0
C              DO JJ = 0, NCS-1
C                JJ1 = JJ
C                JJ2 = JJ1*(JJ-1)
C                JJ3 = JJ2*(JJ-2)
C                JJ4 = JJ3*(JJ-3)
C                Pd0 = Pd0 + CS(JJ)*(dist**JJ)
C                Pd1 = Pd1 + JJ1*CS(JJ)*(dist**(JJ-1))/dls1
C                Pd2 = Pd2 + JJ2*CS(JJ)*(dist**(JJ-2))/dls2
C                Pd3 = Pd3 + JJ3*CS(JJ)*(dist**(JJ-3))/dls3
C                Pd4 = Pd4 + JJ4*CS(JJ)*(dist**(JJ-4))/dls4
C              ENDDO
C
C              IF (Pd0 .GT. 40.D0) THEN           ! Avoid IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
C                G0S = 0.D0                       ! for generating a number < 1.D-38
C                G2S = 0.D0                       ! Second-order derivative for the FF
C                G4S = 0.D0                       ! Fourth-order derivative
C              ELSE IF (Pd0 .LT. -40.D0) THEN
C                G0S = 1.D0
C                G2S = 0.D0
C                G4S = 0.D0
C              ELSE
C                EXPd = EXP(Pd0)
C                G0S  = 1.D0/(1.D0 + EXPd)
C                P1G1 = EXPd*G0S                  ! EXPd*G0S
C                P1G2 = P1G1*G0S                  ! EXPd*(G0S**2)
C                P2G2 = P1G2*EXPd                 ! (EXPd**2)*(G0S**2)
C                P3G3 = P1G1*P2G2                 ! (EXPd**3)*(G0S**3)
C
C                G2S = P1G2*((-1.D0+2.D0*P1G1)*(Pd1**2) - Pd2)
C                G4S = P1G2*((-1.D0+14.D0*P1G1-36.D0*P2G2+24.D0*P3G3)
C     > *(Pd1**4) + (-6.D0+36.D0*P1G1-36.D0*P2G2)*(Pd1**2)*Pd2
C     > + (-4.D0+8.D0*P1G1)*Pd1*Pd3 + (-3.D0+6.D0*P1G1)*(Pd2**2)
C     > - 2.D0*Pd4)
C              ENDIF
C            ELSE
C              G0S = 1.D0
C              G2S = 0.D0
C              G4S = 0.D0
C            ENDIF
C
C            G0(J) = G0E + G0S - 1.D0
C            G2(J) = G2E + G2S
C            G4(J) = G4E + G4S
C          ENDDO
C
        ENDIF
      ENDIF

C      JJ3 = NSTEP/10
C      IF (NSTEP .GT. JJ3*10) JJ3 = JJ3+1
C      DO J = 1, JJ3
C        JJ1 = 1 + (J-1)*10
C        JJ2 = J*10
C        IF (JJ2 .GT. NSTEP) JJ2 = NSTEP
C        WRITE(*,*) 'G0  ', (G0(JJ), JJ = JJ1, JJ2), '  G0'
C        WRITE(*,*) 'G2  ', (G2(JJ), JJ = JJ1, JJ2), '  G2'
C        WRITE(*,*) 'G4  ', (G4(JJ), JJ = JJ1, JJ2), '  G4'
C      ENDDO

C      XLT = XL + XLE + XLS
C      IF(NRES.GT.0) THEN
C        WRITE(NRES,140) PAS, PAS_L, PAS_E, PAS_S, XLT
C      ENDIF

      XLT = XL + XLE + XLS
      IF(NRES.GT.0) THEN
        WRITE(NRES,141) PAS, PAS_L
      ENDIF

  99  FORMAT(/, 5X, 'The  lowest-order nonzero multipole is: ', I0,
     > /, 5X, 'The highest-order nonzero multipole is: ', I0)
 100  FORMAT(/, 5X, ' -----  QUADRUPOLE  : ', 1P,
     > /, 15X,' Length  of  element  = ',G16.8,'  cm',
     > /, 15X,' Bore  radius      RO = ',G13.5,'  cm',
     > /, 15X, 'B-QUADRUPOLE  =', 1P, E15.7, 1X, 'kG',
     > /, 15X, 'The strength of the magnetic quadrupole b2 = '
     >  , E15.7, ' kG/cm')
 101  FORMAT(/, 5X, ' -----  MULTIPOLE   : ', 1P,
     > /, 15X,' Length  of  element  = ',G16.8,'  cm',
     > /, 15X,' Bore  radius      RO = ',G13.5,'  cm',
     > /, 15X, 'B-DIPOLE      =', 1P, E15.7, 1X, 'kG',
     > /, 15X, 'B-QUADRUPOLE  =', 1P, E15.7, 1X, 'kG',
     > /, 15X, 'B-SEXTUPOLE   =', 1P, E15.7, 1X, 'kG',
     > /, 15X, 'B-OCTUPOLE    =', 1P, E15.7, 1X, 'kG',
     > /, 15X, 'B-DECAPOLE    =', 1P, E15.7, 1X, 'kG',
     > /, 15X, 'B-DODECAPOLE  =', 1P, E15.7, 1X, 'kG',
     > /, 15X, 'B-14-POLE     =', 1P, E15.7, 1X, 'kG',
     > /, 15X, 'B-16-POLE     =', 1P, E15.7, 1X, 'kG',
     > /, 15X, 'B-18-POLE     =', 1P, E15.7, 1X, 'kG',
     > /, 15X, 'B-20-POLE     =', 1P, E15.7, 1X, 'kG')

 120  FORMAT(/, 15X, 'Entrance/exit field models are sharp edge',
     >  /, 15X, 'FINTE, FINTS, gap : ', 1P, 3(1X, E12.4))
 121  FORMAT(/, 15X, 'Entrance field model is sharp edge',
     >  /, 15X, 'FINTE, gap : ', 1P, 2(1X, E12.4))
 122  FORMAT(/, 15X, 'Exit field model is sharp edge',
     >  /, 15X, 'FINTS, gap : ', 1P, 2(1X, E12.4))
 128  FORMAT(/,15X,' Entrance  face  ')
 129  FORMAT(/,15X,' Exit  face  ')
 130  FORMAT(20X, ' with  fringe  field :',
     > /, 20X, ' DX  = ', F7.3, '  CM ',
     > /, 20X, ' LAMBDA-QUADRUPOLE =', F7.3, '  CM',
     > /, 20X, I1, ' COEFFICIENTS :', 6F9.5)
 131  FORMAT(20X, 'NOTE: ONLY the sharp-edge model was implemented',
     > ' for the symplectic integrator.' /)
 140  FORMAT(
     >  /, 20X, 'The input  integration step :', 1P, G15.8, ' cm',
     >  /, 20X, 'The ACTUAL integration step :', 1P, G15.8, ' cm',
     >  /, 20X, 'The Entrnc integration step :', 1P, G15.8, ' cm',
     >  /, 20X, 'The Exit   integration step :', 1P, G15.8, ' cm',
     >  /, 20X, 'The total integration length:', 1P, G15.8, ' cm', /)
 141  FORMAT(
     >  /, 20X, 'The input  integration step :', 1P, G15.8, ' cm',
     >  /, 20X, 'The ACTUAL integration step :', 1P, G15.8, ' cm', /)

C---------Compute some constants for the sympletic integration----------
      IF(Q .EQ. 0.D0) Q = 1.0D0
      PREF = BORO*CL9*Q ! the reference momentum; BORO in kG.cm
                        ! CL9 = 0.29979246; Q=1 for proton;
      EREF = SQRT(PREF*PREF + AM*AM) ! the total energy
      BETA0 = PREF / EREF
      GAMMA0 = EREF / AM
      b0g0 = 1.D0/(BETA0*GAMMA0)
      inv_beta0 = 1.D0/BETA0                ! inverse of beta0
      mass = AM / PREF                      ! the scaled mass
      step_L = PAS_L/l0                     ! the scaled stepsize inside
      radius = RO/l0                        ! the scaled radius

      IF (MQI .EQ. 1) THEN                  ! Quadrupole
        eb2 = ((Bfield(2)*RO)/BORO)/(radius**2)   !
        eb(2) = step_L*eb2                  ! Reference momentum P0 = BORO*Q
                                            ! Step_size * eb2/P0
C        IF (NSTEP .GT. NSTEP_L) THEN        ! Presence of fringe field
C          NSTEP_EL = NSTEP_E + NSTEP_L
C          DO J = 1, NSTEP
C            IF (J .LE. NSTEP_E) THEN
C              step = step_E
C            ELSE IF (J .GT. NSTEP_EL) THEN
C              step = step_S
C            ELSE
C              step = step_L
C            ENDIF
C            eb2G0(J) = step*eb2*G0(J)
C            eb2G2(J) = step*eb2*G2(J)/3.D0
C
C            eb2G4(J) = step*eb2*G4(J)/64.D0
C            eb2G4(J) = 0.D0
C          ENDDO

C          JJ3 = NSTEP/10
C          IF (NSTEP .GT. JJ3*10) JJ3 = JJ3+1
C          DO J = 1, JJ3
C            JJ1 = 1 + (J-1)*10
C            JJ2 = J*10
C            IF (JJ2 .GT. NSTEP) JJ2 = NSTEP
C            WRITE(*,*) 'G0  ', (G0(JJ), JJ = JJ1, JJ2), '  G0'
C            WRITE(*,*) 'G2  ', (G2(JJ), JJ = JJ1, JJ2), '  G2'
C          ENDDO

C          IF (newIntegQ .EQ. 1) THEN
C            CALL CHAREF(.FALSE., -XLE, 0.D0, 0.D0)
C---------Option 1: Drift-kick-drift motion integrator for quadrupole---
C            CALL DKD_INTEGR_FF(NSTEP_E, NSTEP_L, NSTEP_S, step_E,
C     >  step_L, step_S, XLE, XLS, b0g0, eb, eb2G0, eb2G2, eb2G4,
C     >  mass, IMAX, F, l0, w0, MQI, HOMP, XLtrc, inv_beta0, G0)
C            CALL CHAREF(.TRUE., -XLS, 0.D0, 0.D0)
C          ELSE
C            CALL ENDJOB ('With fringe field, newIntegQ must be 1', -99)
C          ENDIF

C        ELSE                                ! No fringe field
        IF (newIntegQ .EQ. 1) THEN           ! Only the sharp-edge model
C---------Option 1: Drift-kick-drift motion integrator for quadrupole---
          CALL DKD_INTEGR(NSTEP_L, step_L, b0g0, eb, mass,
     >        IMAX, F, l0, w0, MQI, HOMP)
        ELSE IF (newIntegQ .EQ. 2) THEN
C---------Option 2: Matrix-kick-Matrix motion integrator for quadrupole-
          CALL MKM_INTEGR(NSTEP_L, step_L, b0g0, eb2, mass,
     >        IMAX, F, l0, w0)
        ELSE
          WRITE (*,'(/, 10X, A, I0, A, /)') 'newIntegQ = ',
     >  newIntegQ, ', NOT implemented for quadrupole, stop the job.'
          CALL ENDJOB(' This integrator NOT implemented', -99)
        ENDIF
C        ENDIF

      ELSE IF (MQI .EQ. 2) THEN             ! Multipole
        ROBR = RO/BORO
        DO J = 1, 10
          eb(J) = step_L*(Bfield(J)*ROBR)/(radius**J)    ! step_size * eb_m/P0,
        ENDDO                                            ! where P0 = BORO*Q

        IF (newIntegQ .EQ. 1) THEN
C---------Option 1: Drift-kick-drift motion integrator for multiple-----
          CALL DKD_INTEGR(NSTEP_L, step_L, b0g0, eb, mass,
     >      IMAX, F, l0, w0, MQI, HOMP)
        ELSE
          WRITE (*,'(/, 10X, A, I0, A, /)') 'newIntegQ = ',
     > newIntegQ, ', NOT implemented for multipole, stop the job.'
          CALL ENDJOB(' This integrator has NOT been implemented', -99)
        END IF
      ELSE
        CALL ENDJOB(' Invalid MQI, SHOULD NEVER ENTER HERE', -99)
      ENDIF

      DO IT = 1, IMAX
        IF(NRES.GT.0) THEN
          IF(IT.EQ.1) WRITE(NRES,199)
C----------NOTE: P and T in unit of mrad (as in F(J,I))-----------------
          WRITE(NRES,200) IEX(IT), (FO(J,IT), J=1,5),
     >      XLtrc(IT), (F(J,IT), J=2,6), F(7,IT)*1.D3, IT
        ENDIF
      ENDDO

      CALL SCUMW(XL)
      CALL SCUMR(DUMMY, SCUM, TCUM)
      IF(NRES .GT. 0) THEN
        WRITE(NRES,201) KP, XCE, YCE, ALE, SCUM*UNIT(5), TCUM
      ENDIF

C      ppii = 13591409.D0 - 120.D0*
C     >  (13591409.D0 + 545140134.D0)/(640320.D0**3)
C      ppii = ppii + 12*11*10*9*8*7*(13591409.D0 + 545140134.D0*2)/
C     > (8.D0*(640320.D0**6))
C      ppii = (12.D0/(640320**1.5D0)) * ppii
C      ppii = 1.D0/ppii
C      write (*,*) 'ppii = ', ppii

 199  FORMAT(2X, '  KPOS  DP         Y(cm)   T(mrad)     ',
     >  'Z(cm)   P(mrad)   |', 2X, ' X(cm)         Y(cm)        ',
     >  'T(mrdd)      Z(cm)        P(mrad)     S(cm)         Time(ns)',
     >   4X, '  I')
 200  FORMAT(2X,'A',2X,I3,F8.4,4F10.3, '   |', F12.6, 6F13.6, 2X, I5)
 201  FORMAT(/, 5X, 'KPOS =  ', I0,
     >  '.  Change  of  frame  at  exit  of  element.', /, 10X,
     >  'X =', 1P, G12.4, ' CM   Y =', G12.4, 1P,
     >  ' cm,  tilt  angle =', G14.6, ' RAD', /, /, /, 1P,
     >  ' Cumulative length of optical axis = ', 1P, G17.9, ' m ;  ',
     >  'Time  (for ref. rigidity & particle) = ', 1P, G14.6, ' s ')

      RETURN
      END
