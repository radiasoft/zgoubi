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
      SUBROUTINE DKD_INTEGR_FF(NSTEP_E, NSTEP_L, NSTEP_S, step_E,
     >  step_L, step_S, XLE, XLS, b0g0, eb, eb2G0, eb2G2, eb2G4,
     >  mass, IMAX, data, l0, w0, MQI, HOMP, XLtrc, inv_beta0)
C     --------------------------------------------------------------------
C     Drift(L/2)Kick(L)Drift(L/2) motion integrator for
C           both magnetic quadrupole and multipole
C     ALL the variables are scaled in this subroutine
C     --------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE "MAXTRA.H"      ! PARAMETER (MXJ=7)
      INCLUDE "MAXCOO.H"      ! PARAMETER (MXT=10000)
      INCLUDE "MXSTEP.H"      ! PARAMETER (MXSTEP=1000000)

      INTEGER, INTENT(IN) :: NSTEP_E, NSTEP_L, NSTEP_S
      INTEGER, INTENT(IN) :: IMAX, MQI, HOMP
      DOUBLE PRECISION, INTENT(IN) :: step_E, step_L, step_S, XLE
      DOUBLE PRECISION, INTENT(IN) :: b0g0, mass, l0, w0, XLS, inv_beta0
      DOUBLE PRECISION, INTENT(IN) :: eb(10),eb2G0(MXSTEP),eb2G2(MXSTEP)
      DOUBLE PRECISION, INTENT(IN) :: eb2G4(MXSTEP)
      DOUBLE PRECISION, INTENT(IN OUT) :: data(MXJ,MXT)
      DOUBLE PRECISION, INTENT(OUT) :: XLtrc(MXT)
      DOUBLE PRECISION step, halfStep_E, halfStep_L, halfStep_S
      DOUBLE PRECISION Ptot, theta, phi, X, Y, PX, PY, T, PT, S
      DOUBLE PRECISION P2, Ps, Ps2, Pxy2, Pxz, delT, alpha
      DOUBLE PRECISION XC, ZC, DX, DL, DS, DT, XLT
      DOUBLE PRECISION alpha_ff1, alpha_ff2, alpha_ff3, alpha_ff4
      DOUBLE PRECISION Pdelta, Pdelta3, PTbeta, X3, Y3
      DOUBLE PRECISION X0, Y0, PX0, PY0, cospi4, sinpi4
      INTEGER IT, ISTEP, NSTEP, NSTEP_EL

C--------Zgoubi uses (Y, T, Z, P, SAR, TAR) coordinates-----------------
C      DP=F(1,I)
C      Y =F(2,I)
C      T =F(3,I)*0.001D0
C      Z =F(4,I)
C      P =F(5,I)*0.001D0
C      SAR=F(6,I)
C      TAR=F(7,I)*1.D5
C-----------------------------------------------------------------------

C-----------Constants for all the particles-----------------------------
      halfStep_L = 0.5D0*step_L       ! 0.5L/l0, scaled length
      halfStep_E = 0.5D0*step_E
      halfStep_S = 0.5D0*step_S
      NSTEP_EL = NSTEP_E + NSTEP_L
      NSTEP = NSTEP_EL + NSTEP_S
      alpha_ff1 = eb(2)/6.D0
      alpha_ff2 = eb(2)/12.D0
      alpha_ff3 = eb(2)/2.D0
      alpha_ff4 = eb(2)/4.D0
      cospi4 = SQRT(0.5D0)
      sinpi4 = cospi4

      DO CONCURRENT (IT = 1:IMAX)     ! Loop over all the particles
        Ptot = data(1,IT)             ! the scaled total momentum
        PT = SQRT(Ptot**2 + mass**2)  ! the scaled total energy
        P2 = PT*PT - b0g0*b0g0        ! P^2 = PT^2 - (beta0*gamma0)^(-2)
        Pdelta = SQRT(P2)
        Pdelta3 = Pdelta**3
        PTbeta = PT + inv_beta0

        phi   = data(5,IT)*0.001D0    ! angle phi (in rad)
        theta = data(3,IT)*0.001D0    ! angle theta (in rad)

C---------Convert the coordinates from Zgobi notation to DTA notation,
C---------Frame change by -XLE (in DTA notation, Z <----> X of Zgoubi
C                                                X <----> Y of Zgoubi
C                                                Y <----> Z of Zgoubi
        X = data(2,IT)
        Y = data(4,IT)
        S = data(6,IT)
        T = data(7,IT)*1.D-6

        ZC = -XLE                     ! ZC <----> XC of Zgoubi
        XC = 0.D0                     ! XC <----> YC of Zgoubi
        alpha = 0.D0                  ! alpha <-> A  of Agoubi

        DX = ((X-XC)*COS(theta) + ZC*SIN(theta))/COS(theta-alpha) - X
        DL = SQRT(DX*DX + ZC*ZC)
        DL = SIGN(DL, ZC)
        DS = DL/COS(phi)
        DT = DS / ((Ptot/PT)*(w0*l0)) ! l0*w0 = speed of light in cm/s

        X = X + DX
        Y = Y + DL*TAN(phi)
        S = S + DS
        T = T + DT

C---------then scaled by l0, w0, and P0 = PREF--------------------------
        X = X/l0
        Y = Y/l0
        S = S/l0
        T = T*(-w0)

        PX= Ptot*cos(phi)*sin(theta)      ! Ptot is already scaled by PREF
        PY= Ptot*sin(phi)
        XLT = 0.D0
C---------Start the symplectic dfift-kick-drift integrator--------------
        ISTEP = 1
        DO 999 WHILE (ISTEP .LE. NSTEP)
C----------DRIFT: half the stepsize for the first drift-----------------
          IF ((ISTEP .EQ. 1) .OR. (ISTEP .EQ. NSTEP_E+1) .OR.
     >     (ISTEP .EQ. NSTEP_EL+1)) THEN

            IF (ISTEP .EQ. 1) THEN
              step = halfStep_E
              IF (NSTEP_E .EQ. 0) step = halfStep_L
            ENDIF

            IF (ISTEP .EQ. NSTEP_E +1) step = halfStep_L
            IF (ISTEP .EQ. NSTEP_EL+1) step = halfStep_S

            Pxy2 = PX*PX + PY*PY
            Ps2  = P2 - Pxy2
            Ps   = SQRT(Ps2)
            delT = step/Ps                ! the scaled drift time

            X = X + delT*PX
            Y = Y + delT*PY
            S = S + step*SQRT(1+Pxy2/Ps2)
            T = T - delT*PT
            XLT = XLT + step
          ENDIF

C----------KICK---------------------------------------------------------
          IF (MQI .EQ. 1) THEN            ! Quadrupole

            PX = PX - eb2G0(ISTEP)*X      ! Soft fringe-field is treated as
            PY = PY + eb2G0(ISTEP)*Y      ! the zeroth-order perturbation

C            PX = PX - eb2G0(ISTEP)*X + eb2G2(ISTEP)*(X**3)
C     > - eb2G4(ISTEP)*((X**5)/2.D0+(X**3)*(Y**2)/3.D0-X*(Y**4)/6.D0)
C             PY = PY + eb2G0(ISTEP)*Y - eb2G2(ISTEP)*(Y**3)
C     > + eb2G4(ISTEP)*((Y**5)/2.D0+(Y**3)*(X**2)/3.D0-Y*(X**4)/6.D0)

          ELSE IF (MQI .EQ. 2) THEN       ! Multiple
            PX = PX - eb(1) - eb(2)*X - eb(3)*(X*X-Y*Y)
            PY = PY + eb(2)*Y + eb(3)*2.D0*X*Y

            IF(HOMP .GT. 3) THEN          ! Rarely occur when highest-order
                                          ! nonezero multipole > 3
              IF (eb(4) .NE. 0.D0) THEN   ! Contributions from octupole
                PX = PX - eb(4)*(X*X*X - 3.D0*X*Y*Y)
                PY = PY + eb(4)*(3.D0*X*X*Y - Y*Y*Y)
              ENDIF

              IF (eb(5) .NE. 0.D0) THEN   ! Contributions from decapole
                PX = PX - eb(5)*(X**4 - 6.D0*X*X*Y*Y + Y**4)
                PY = PY + eb(5)*(4.D0*X*X*X*Y - 4.D0*X*Y*Y*Y)
              ENDIF

              IF (eb(6) .NE. 0.D0) THEN   ! Contributions from dodacapole
                PX = PX - eb(6)*(X**5 -10.D0*X*X*X*Y*Y +5.D0*X*(Y**4))
                PY = PY + eb(6)*(5.D0*(X**4)*Y -10.D0*X*X*Y*Y*Y +Y**5)
              ENDIF

              IF (eb(7) .NE. 0.D0) THEN   ! Contributions from 14-pole
                PX = PX-eb(7)*(X**6-15.D0*(X**4)*Y*Y+15.D0*X*X*(Y**4))
                PY = PY+eb(7)*(6.D0*(X**5)*Y - 20.D0*X*X*X*Y*Y*Y
     >                                       + 6.D0*X*Y**5)
              ENDIF

              IF (eb(8) .NE. 0.D0) THEN   ! Contributions from 16-pole
                PX = PX - eb(8)*(X**7 - 21.D0*(X**5)*Y*Y
     >                  + 35.D0*(X**3)*(Y**4) - 7.D0*X*(Y**6))
                PY = PY + eb(8)*(7.D0*(X**6)*Y - 35.D0*(X**4)*Y*Y*Y
     >                  + 21.D0*X*X*(Y**5) - Y**7)
              ENDIF

              IF (eb(9) .NE. 0.D0) THEN   ! Contributions from 18-pole
                PX = PX - eb(9)*(X**8 - 28.D0*(X**6)*Y*Y
     >             + 70.D0*(X**4)*(Y**4) - 28.D0*X*X*(Y**6) + Y**8)
                PY = PY + eb(9)*(8.D0*(X**7)*Y - 56.D0*(X**5)*Y*Y*Y
     >                  + 56.D0*X*X*X*(Y**5) - 8.D0*X*(Y**7))
              ENDIF

              IF (eb(10).NE. 0.D0) THEN   ! Contributions from 20-pole
                PX = PX - eb(10)*(X**9 - 36.D0*(X**7)*Y*Y
     >                  + 126.D0*(X**5)*(Y**4) - 84.D0*X*X*X*(Y**6)
     >                  + 9.D0*Y**9)
                PY = PY + eb(10)*(9.D0*(X**8)*Y - 84.D0*(X**6)*Y*Y*Y
     >             + 126.D0*(X**4)*(Y**5) - 36.D0*X*X*(Y**7) + Y**9)
              ENDIF
            ENDIF
          ENDIF

C----------DRIFT: the whole step size except for the last step----------
          IF (ISTEP .LE. NSTEP_E) THEN
            step = step_E
          ELSE IF (ISTEP .GT. NSTEP_EL) THEN
            step = step_S
          ELSE
            step = step_L
          ENDIF

          IF (ISTEP .EQ. NSTEP) THEN  ! Stepsize = L/2 for the last step
            step = halfStep_S
            IF (NSTEP_S .EQ. 0) step = halfStep_L
          ENDIF

          IF (ISTEP .EQ. NSTEP_E)  step = halfStep_E
          IF (ISTEP .EQ. NSTEP_EL) step = halfStep_L

          Pxy2 = PX*PX + PY*PY
          Ps2  = P2 - Pxy2
          Ps   = SQRT(Ps2)
          delT = step/Ps

          X = X + delT*PX
          Y = Y + delT*PY
          S = S + step*SQRT(1+Pxy2/Ps2)
          T = T - delT*PT
          XLT = XLT + step

          ISTEP = ISTEP + 1

C--------KICK due to the hard-edge fringe field model-------------------
          IF (MQI .EQ. 1) THEN                  ! Quadrupole
            IF (ISTEP .EQ. (NSTEP_E+1)) THEN    ! Entrance of the fringe field, the hard-edge model
              X0  = X
              Y0  = Y
              PX0 = PX
              PY0 = PY
              X  =  X0*cospi4 -  Y0*sinpi4
              Y  =  X0*sinpi4 +  Y0*cospi4
              PX = PX0*cospi4 - PY0*sinpi4
              PY = PX0*sinpi4 + PY0*cospi4

              X3 = X*X*X
              PX = PX + alpha_ff4*X*X*PY/Pdelta
              Y  = Y  - alpha_ff2*X3/Pdelta
              T  = T  + alpha_ff2*X3*PY*PTbeta/Pdelta3

              Y3 = Y*Y*Y
              PY = PY + alpha_ff3*Y*Y*PX/Pdelta
              X  = X  - alpha_ff1*Y3/Pdelta
              T  = T  + alpha_ff1*Y3*PX*PTbeta/Pdelta3

              X3 = X*X*X
              PX = PX + alpha_ff4*X*X*PY/Pdelta
              Y  = Y  - alpha_ff2*X3/Pdelta
              T  = T  + alpha_ff2*X3*PY*PTbeta/Pdelta3

              X0  = X
              Y0  = Y
              PX0 = PX
              PY0 = PY
              X  =  X0*cospi4 +  Y0*sinpi4
              Y  = -X0*sinpi4 +  Y0*cospi4
              PX = PX0*cospi4 + PY0*sinpi4
              PY =-PX0*sinpi4 + PY0*cospi4
            ENDIF

            IF (ISTEP .EQ. (NSTEP_EL+1)) THEN   ! Exit of the fringe field, the hard-edge model
              X0  = X
              Y0  = Y
              PX0 = PX
              PY0 = PY
              X  =  X0*cospi4 -  Y0*sinpi4
              Y  =  X0*sinpi4 +  Y0*cospi4
              PX = PX0*cospi4 - PY0*sinpi4
              PY = PX0*sinpi4 + PY0*cospi4

              X3 = X*X*X
              PX = PX - alpha_ff4*X*X*PY/Pdelta
              Y  = Y  + alpha_ff2*X3/Pdelta
              T  = T  - alpha_ff2*X3*PY*PTbeta/Pdelta3

              Y3 = Y*Y*Y
              PY = PY - alpha_ff3*Y*Y*PX/Pdelta
              X  = X  + alpha_ff1*Y3/Pdelta
              T  = T  - alpha_ff1*Y3*PX*PTbeta/Pdelta3

              X3 = X*X*X
              PX = PX - alpha_ff4*X*X*PY/Pdelta
              Y  = Y  + alpha_ff2*X3/Pdelta
              T  = T  - alpha_ff2*X3*PY*PTbeta/Pdelta3

              X0  = X
              Y0  = Y
              PX0 = PX
              PY0 = PY
              X  =  X0*cospi4 +  Y0*sinpi4
              Y  = -X0*sinpi4 +  Y0*cospi4
              PX = PX0*cospi4 + PY0*sinpi4
              PY =-PX0*sinpi4 + PY0*cospi4
             ENDIF
           ENDIF

 999    CONTINUE

C----------Revert coordinates from DTA to Zgoubi------------------------
        Pxz = SQRT(Ptot*Ptot - PY*PY) ! the projected momentum on the xz-plane
        theta = ASIN(PX/Pxz)          ! theta = arcsin(px/p_y)
        phi = ASIN(PY/Ptot)           ! phi = arcsin(py/P)
C----------Frame Change by -XLS due to exit fringe field----------------
        ZC = -XLS/l0                  ! Here X, Y, S, T are all scaled
        XC = 0.D0/l0                  ! So we need to scale ZC, XC
        alpha = 0.D0

        DX = ((X-XC)*COS(theta) + ZC*SIN(theta))/COS(theta-alpha) - X
        DL = SQRT(DX*DX + ZC*ZC)
        DL = SIGN(DL, ZC)
        DS = DL/COS(phi)
        DT = DS / ((Ptot/PT)*w0)   ! Get the time in unit of s

        X = X + DX
        Y = Y + DL*TAN(phi)
        S = S + DS
        T = T + DT*(-w0)                     ! T is scaled, but DT is NOT
C----------Scale them back----------------------------------------------
        data(2,IT) = X*l0
        data(3,IT) = theta*1000.D0           ! theta in mrad in F-matrix
        data(4,IT) = Y*l0
        data(5,IT) = phi*1000.D0             ! phi in mrad in F-matrix
        data(6,IT) = S*l0
        data(7,IT) = -T*1.D6/w0              ! time in us in F-matrix
        XLtrc(IT)  = XLT*l0                  ! Make sure the ray-tracing
                                             ! length = XL + XLE + XLS
      ENDDO

      RETURN
      END
