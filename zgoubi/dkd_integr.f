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
      SUBROUTINE DKD_INTEGR(MAX_STEP, stepsize, b0g0, eb2, mass,
     >      IMAX, data, X_axis)
C     --------------------------------------------------------------------
C     Drift(L/2)Kick(L)Drift(L/2) motion integrator for magnetic quadrupole
C     ALL the variables are scaled
C     --------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE "MAXTRA.H"      ! PARAMETER (MXJ=7)
      INCLUDE "MAXCOO.H"      ! PARAMETER (MXT=10000)
      INCLUDE "SCALE_symp.H"  ! PARAMETER (l0, w0)

      INTEGER, INTENT(IN) :: MAX_STEP, IMAX
      DOUBLE PRECISION, INTENT(IN) :: stepsize, b0g0, eb2, mass
      DOUBLE PRECISION, INTENT(IN OUT) :: data(MXJ,MXT)
      DOUBLE PRECISION, INTENT(OUT) :: X_axis(MXT)
      DOUBLE PRECISION halfStep, KL, Ptot, step0
      DOUBLE PRECISION theta, phi, X, Y, Z, PX, PY, T, PT, S
      DOUBLE PRECISION P2, Ps, Ps2, Pxy2, Pxz, delT
      INTEGER NUM_STEP, IT

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
      halfStep = 0.5D0*stepsize                ! 0.5L/L0
      KL = eb2*stepsize

      DO CONCURRENT (IT = 1:IMAX)     ! Loop over all the particles
        Ptot = data(1,IT)             ! the scaled total momentum
        PT = SQRT(Ptot**2 + mass**2)  ! the scaled total energy
        phi   = data(5,IT)*0.001D0    ! angle phi
        theta = data(3,IT)*0.001D0    ! angle theta

C---------Convert the coordinates from Zgobi notation to DTA notation,
C---------then scaled by l0, w0, and P0 = PREF--------------------------
        Z = data(6,IT)*cos(phi)*cos(theta)/l0
        X = data(2,IT)/l0
        Y = data(4,IT)/l0
        PX= Ptot*cos(phi)*sin(theta)  ! Ptot is already scaled by PREF
        PY= Ptot*sin(phi)
        T = data(7,IT)*(-w0)          ! time = F(7,I), but TAR = F(7,I)*1.D5
        S = data(6,IT)/l0             ! displacement = F(6,I)

C---------Compute constants for a particle------------------------------
        P2 = PT*PT - b0g0*b0g0

C---------Start the symplectic dfift-kick-drift integrator--------------
        NUM_STEP = 1
        DO 999 WHILE (NUM_STEP .LE. MAX_STEP)
C----------DRIFT: half the stepsize for the first drift-----------------
          IF (NUM_STEP .EQ. 1) THEN
            Pxy2 = PX*PX + PY*PY
            Ps2  = P2 - Pxy2
            Ps   = SQRT(Ps2)
            delT = halfStep/Ps        ! the scaled drift time

            X = X + delT*PX
            Y = Y + delT*PY
            T = T - delT*PT
            Z = Z + halfStep*SQRT(1-Pxy2/Ps2)
            S = S + halfStep
          ENDIF

C----------KICK---------------------------------------------------------
          PX = PX - KL*X
          PY = PY + KL*Y

C----------DRIFT: the whole step size except for the last step----------
          step0 = stepsize
          IF (NUM_STEP .EQ. MAX_STEP) step0 = halfStep ! Stepsize = L/2 for the last step

          Pxy2 = PX*PX + PY*PY
          Ps2  = P2 - Pxy2
          Ps   = SQRT(Ps2)
          delT = step0/Ps

          X = X + delT*PX
          Y = Y + delT*PY
          T = T - delT*PT
          Z = Z + step0*SQRT(1-Pxy2/Ps2)
          S = S + step0

          NUM_STEP = NUM_STEP + 1
 999    CONTINUE

C----------Revert coordinates from DTA to Zgoubi------------------------
C----------And then scale them back-------------------------------------
        Pxz = SQRT(Ptot*Ptot - PY*PY) ! the projected momentum on the xz-plane
        theta = ASIN(PX/Pxz)  ! theta = arcsin(px/p_y)
        phi = ASIN(PY/Ptot)   ! phi = arcsin(py/P)

        X_axis(IT) = Z*l0
        data(2,IT) = X*l0
        data(3,IT) = theta*1000.D0
        data(4,IT) = Y*l0
        data(5,IT) = phi*1000.D0
        data(6,IT) = S*l0
        data(7,IT) = -T/w0

      ENDDO

      RETURN
      END
