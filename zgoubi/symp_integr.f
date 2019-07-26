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
C     (3) The drift-kick-drift symplectic particle motion integrator
C         for a general magnetic multipole (MQchoice = 2, newIntegQ = 1)
C-----------------------------------------------------------------------
      INCLUDE "MAXTRA.H"   ! PARAMETER (MXJ=7)
      INCLUDE "MAXCOO.H"   ! PARAMETER (MXT=10000)
      INCLUDE "MXFS.H"     ! PARAMETRE (MXF=45, MXS=2000, MLF=9, MXP=11)
      INCLUDE "MXLD.H"     ! PARAMETER (MXL=15000 , MXD=1400)
      INCLUDE "MXSCL.H"    ! PARAMETER (MXSCL=10)
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

      DOUBLE PRECISION XL, XE, RO, GAP, ROBR, XS, YS, ZR, ZS, YR
      DOUBLE PRECISION DLE, XLS, DLS, DL0, TCUM, SCUM, l0, w0, eb2
      DOUBLE PRECISION FINTE, FINTS, PREF, EREF, BETA0, GAMMA0
      DOUBLE PRECISION b0g0, PAS0, step, radius, mass, DUMMY
      DOUBLE PRECISION Bfield(10), RSA(10), eb(10)
      DOUBLE PRECISION FFE(10), FFCE(0:5), FFS(10), FFCS(0:5)
      INTEGER KUASEX, NRES, newIntegQ, NSTEP, J,IT,NCE,NCS, MQI, HOMP

      PARAMETER (l0 = 100.D0)
      PARAMETER (w0 = 299792458.D0)

C-------- Get the parameters for magnetic quadrupole (MQI .EQ. 1) ------
C---------or magnetic multipole (MQI .EQ. 2) ---------------------------
      IF (MQI .EQ. 1) THEN
        XL = A(NOEL,10)
        RO = A(NOEL,11)
        Bfield(2) = A(NOEL,12)

        XE = A(NOEL,20)                      ! Integration zone: entrance
        DLE= A(NOEL,21)                      ! Dipole fringe field lambda_E
        NCE = NINT(A(NOEL,30))               ! NCE = unused number
        XLS = A(NOEL,40)                     ! Integration zone: exit
        DLS = A(NOEL,41)                     ! Dipole fringe field lambda_S
        NCS = NINT(A(NOEL,50))               ! NCS = unused number
        DO J = 0, 5
          FFCE(J) = A(NOEL,31+J)             ! Entrance fringe field coefficients
          FFCS(J) = A(NOEL,51+J)             ! Exit fringe field coefficients
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
            FFE(J) = A(NOEL,14+J)            ! Quadrupole fringe field = lambda_E * FFE(J != 1)
            FFS(J) = A(NOEL,32+J)            ! Quadrupole fringe field = lambda_S * FFS(J != 1)
          ENDIF
          IF (J .LE. 6) THEN
            FFCE(J-1) = A(NOEL,25+J)         ! Entrance fringe field coefficients
            FFCS(J-1) = A(NOEL,43+J)         ! Exit fringe field coefficients
          ENDIF
        ENDDO

        XE = A(NOEL,14)                      ! Integration zone: entrance
        DLE= A(NOEL,15)                      ! Dipole fringe field lambda_E
        NCE = NINT(A(NOEL,25))               ! NCE = unused number
        XLS = A(NOEL,32)                     ! Integration zone: exit
        DLS = A(NOEL,33)                     ! Dipole fringe field lambda_S
        NCS = NINT(A(NOEL,43))               ! NCS = unused number

        PAS = A(NOEL,60)                     ! Integration step
        KP = NINT(A(NOEL,63))                ! In the Zgoubi manu, KPOS is located at # 61
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
          WRITE(NRES, 101) XL, RO, (Bfield(J), J=1,10)
        ENDIF

        KUASEX = 0
        HOMP = 0
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

      DL0 = DLE + DLS
      IF(DL0 .EQ. 0.D0) THEN
C-------- Sharp edge at entrance and exit---------------------
        FINTE = XE
        XE = 0.D0
        FINTS = XLS
        XLS = 0.D0
        IF(NRES.GT.0) THEN
          WRITE(NRES,102) FINTE, FINTS, GAP
        ENDIF
      ENDIF

      PAS0 = PAS                  ! Integration step from input
      NSTEP = CEILING(XL/PAS0)
      PAS = XL/NSTEP              ! Adjust the step size to fit the length
      IF(NRES.GT.0) THEN
        WRITE(NRES,103) PAS0, PAS
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

 102  FORMAT(/, 15X, 'Entrance/exit field models are sharp edge',
     >  /, 15X, 'FINTE, FINTS, gap : ', 1P, 3(1X, E12.4))
 103  FORMAT(
     >  /, 20X, 'The input  integration step :', 1P, G15.8, ' cm',
     >  /, 20X, 'The ACTUAL integration step :', 1P, G15.8, ' cm', /)

C---------Compute some constants for the sympletic integration----------
      IF(Q .EQ. 0.D0) Q = 1
      PREF = BORO*CL9*Q ! the reference momentum; BORO in kG.cm
                        ! CL9 = 0.29979246; Q=1 for proton;
      EREF = SQRT(PREF*PREF + AM*AM) ! the total energy
      BETA0 = PREF / EREF
      GAMMA0 = EREF / AM
      b0g0 = 1.D0/(BETA0*GAMMA0)
      mass = AM / PREF                      ! the scaled mass
      step = PAS/l0                         ! the scaled stepsize
      radius = RO/l0                        ! the scaled radius

      IF (MQI .EQ. 1) THEN                  ! Quadrupole
        eb2 = ((Bfield(2)*RO)/BORO)/(radius**2)   ! reference momentum P0 = BORO*Q
        eb(2) = step*eb2                          ! step_size * eb2/P0

        IF (newIntegQ .EQ. 1) THEN
C---------Option 1: Drift-kick-drift motion integrator for quadrupole---
          CALL DKD_INTEGR(NSTEP, step, b0g0, eb, mass,
     >                           IMAX, F, l0, w0, MQI, HOMP)
        ELSE IF (newIntegQ .EQ. 2) THEN
C---------Option 2: Matrix-kick-Matrix motion integrator for quadrupole-
          CALL MKM_INTEGR(NSTEP,step,b0g0,eb2,mass,IMAX,F,l0,w0)
        ELSE
          WRITE (*,'(/, 10X, A, I0, A, /)') 'newIntegQ = ',
     >  newIntegQ, ', NOT implemented for quadrupole, stop the job.'
          CALL ENDJOB(' This integrator has NOT been implemented', -99)
        END IF

      ELSE IF (MQI .EQ. 2) THEN             ! Multipole
        ROBR = RO/BORO
        DO J = 1, 10
          eb(J) = step*(Bfield(J)*ROBR)/(radius**J) ! step_size * eb_m/P0,
        ENDDO                                       ! where P0 = BORO*Q

        IF (newIntegQ .EQ. 1) THEN
C---------Option 1: Drift-kick-drift motion integrator for multiple-----
          CALL DKD_INTEGR(NSTEP, step, b0g0, eb, mass,
     >                           IMAX, F, l0, w0, MQI, HOMP)
        ELSE
          WRITE (*,'(/, 10X, A, I0, A, /)') 'newIntegQ = ',
     > newIntegQ, ', NOT implemented for multipole, stop the job.'
          CALL ENDJOB(' This integrator has NOT been implemented', -99)
        END IF
      ELSE
        CALL ENDJOB(' Invalid MQI', -99)
      ENDIF

      DO IT = 1, IMAX
        IF(NRES.GT.0) THEN
          IF(IT.EQ.1) WRITE(NRES,199)
C----------NOTE: P and T in unit of mrad (as in F(J,I))-----------------
          WRITE(NRES,200) IEX(IT), (FO(J,IT), J=1,5),
     >      XL, (F(J,IT), J=2,6), F(7,IT)*1.D3, IT
        ENDIF
      ENDDO

      CALL SCUMW(XL)
      CALL SCUMR(DUMMY, SCUM, TCUM)
      IF(NRES .GT. 0) THEN
        WRITE(NRES,201) KP, XCE, YCE, ALE, SCUM*UNIT(5), TCUM
      ENDIF

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
