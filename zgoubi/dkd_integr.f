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
      SUBROUTINE DKD_INTEGR(MAX_STEP, PAS, b0g0, eb2, AM, PREF, IMAX,
     >      data, X_axis)
C     --------------------------------------------------------------------
C     Drift(L/2)Kick(L)Drift(L/2) motion integrator for magnetic quadrupole
C     ALL the variables are scaled
C     --------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE "MAXTRA.H"      ! PARAMETER (MXJ=7)
      INCLUDE "MAXCOO.H"      ! PARAMETER (MXT=10000)
      INCLUDE "SCALE_symp.H"  ! PARAMETER (l0, w0)

      INTEGER, INTENT(IN) :: MAX_STEP, IMAX
      DOUBLE PRECISION, INTENT(IN) :: PAS, b0g0, eb2, PREF, AM
      DOUBLE PRECISION, INTENT(IN OUT) :: data(MXJ,MXT)
      DOUBLE PRECISION, INTENT(OUT) :: X_axis(MXT)
      DOUBLE PRECISION theta, phi, x, y, z, s, t, px, py, pt
      DOUBLE PRECISION S_scl, Z_scl, X_scl, PX_scl, Y_scl
      DOUBLE PRECISION PY_scl, T_scl, PT_scl, KL, PAS0, hfPAS
      DOUBLE PRECISION P2, Ps, Ps2, Pxy2, Tdft, Pxz, Ptot, Etot
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
      hfPAS = 0.5D0*PAS                ! 0.5L/L0C      DP=F(1,I)
      KL = eb2*PAS

      DO CONCURRENT (IT = 1:IMAX)  !!! Loop over all the particles
        Ptot = PREF*data(1,IT)         ! the total momentum
        Etot = SQRT(Ptot*Ptot + AM*AM) ! the total energy
        phi   = data(5,IT)*0.001D0     ! angle phi
        theta = data(3,IT)*0.001D0     ! angle theta

C---------Convert the coordinate system from Zgobi to DTA---------------
        z = data(6,IT)*cos(phi)*cos(theta)
        x = data(2,IT)
        y = data(4,IT)
        px= Ptot*cos(phi)*sin(theta)
        py= Ptot*sin(phi)
        pt= Etot
        t = data(7,IT)        ! time = F(7,I), but TAR = F(7,I)*1.D5
        s = data(6,IT)        ! displacement = F(6,I)

C---------Scale the coordinated from DTA to dimensionless---------------
        X_scl = x/l0
        Y_scl = y/l0
        Z_scl = z/l0
        PX_scl = px/PREF
        PY_scl = py/PREF
        PT_scl = pt/PREF
        T_scl = -t*w0
        S_scl =  s/l0

C---------Compute constants for a particle------------------------------
        P2 = PT_scl*PT_scl - b0g0*b0g0

C---------Start the symplectic dfift-kick-drift integrator--------------
        NUM_STEP = 1
        DO 999 WHILE (NUM_STEP .LE. MAX_STEP)
C----------DRIFT: half the stepsize for the first drift-----------------
          IF (NUM_STEP .EQ. 1) THEN
            Pxy2 = PX_scl*PX_scl + PY_scl*PY_scl
            Ps2  = P2 - Pxy2
            Ps   = SQRT(Ps2)
            Tdft = hfPAS/Ps                      ! the scaled drift time

            X_scl = X_scl + Tdft*PX_scl
            Y_scl = Y_scl + Tdft*PY_scl
            T_scl = T_scl - Tdft*PT_scl
            Z_scl = Z_scl + hfPAS*SQRT(1-Pxy2/Ps2)
            S_scl = S_scl + hfPAS
          ENDIF

C----------KICK---------------------------------------------------------
          PX_scl = PX_scl - KL*X_scl
          PY_scl = PY_scl + KL*Y_scl

C----------DRIFT: the whole step size except for the last step----------
          PAS0 = PAS
          IF (NUM_STEP .EQ. MAX_STEP) PAS0 = hfPAS ! Stepsize = L/2 for the last step

          Pxy2 = PX_scl*PX_scl + PY_scl*PY_scl
          Ps2  = P2 - Pxy2
          Ps   = SQRT(Ps2)
          Tdft = PAS0/Ps

          X_scl = X_scl + Tdft*PX_scl
          Y_scl = Y_scl + Tdft*PY_scl
          T_scl = T_scl - Tdft*PT_scl
          Z_scl = Z_scl + PAS0*SQRT(1-Pxy2/Ps2)
          S_scl = S_scl + PAS0

          NUM_STEP = NUM_STEP + 1
 999    CONTINUE

C----------Revert scaling coordinates from dimensionless to DTA---------
        x = X_scl*l0
        y = Y_scl*l0
        z = Z_scl*l0
        px= PX_scl*PREF
        py= PY_scl*PREF
        t = -T_scl/w0
        s = S_scl*l0

C----------Revert coordinates from DTA to Zgoubi------------------------
C----------Assuming phi and theta are both quite small------------------
        Pxz = SQRT(Ptot*Ptot - py*py) ! the projected momentum on the xz-plane
        theta = ASIN(px/Pxz)  ! theta = arcsin(px/p_y)
        phi = ASIN(py/Ptot)   ! phi = arcsin(py/P)

        X_axis(IT) = z
        data(2,IT) = x
        data(3,IT) = theta*1000.D0
        data(4,IT) = y
        data(5,IT) = phi*1000.D0
        data(6,IT) = s
        data(7,IT) = t

      ENDDO

      RETURN
      END
