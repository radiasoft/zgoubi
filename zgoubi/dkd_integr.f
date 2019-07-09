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
      SUBROUTINE DKD_INTEGR(MAX_STEP, PAS_scl, hfPAS_scl,
     >      gb0, eb2, PREF, IMAX, AM, data, X_axis)
C     --------------------------------------------------------------------
C     Drift(L/2)Kick(L)Drift(L/2) motion integrator for magnetic quadrupole
C     ALL the variables are scaled
C     --------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE "MAXTRA.H"      ! PARAMETER (MXJ=7)
      INCLUDE "MAXCOO.H"      ! PARAMETER (MXT=10000)
      INCLUDE "SCALE_symp.H"  ! PARAMETER (scl_l0, scl_w0)

      INTEGER, INTENT(IN) :: MAX_STEP, IMAX
      DOUBLE PRECISION, INTENT(IN) :: PAS_scl,hfPAS_scl,gb0,eb2,PREF,AM
      DOUBLE PRECISION, INTENT(IN OUT) :: data(MXJ,MXT)
      DOUBLE PRECISION, INTENT(OUT) :: X_axis(MXT)
      DOUBLE PRECISION theta, phi, x_dta, y_dta, z_dta, s_dta, t_dta
      DOUBLE PRECISION px_dta, py_dta, pt_dta, S_scl, Z_scl
      DOUBLE PRECISION X_scl, PX_scl, Y_scl, PY_scl, T_scl, PT_scl
      DOUBLE PRECISION pt_gb2, ps, portion, deltaZ, DTscl
      DOUBLE PRECISION Pxz, Ptot, Etot
      INTEGER NUM_STEP, IT

C      DP=F(1,I)
C      Y =F(2,I)
C      T =F(3,I)*0.001D0
C      Z =F(4,I)
C      P =F(5,I)*0.001D0
C      SAR=F(6,I)
C      TAR=F(7,I)*1.D5

      DO CONCURRENT (IT = 1:IMAX)  !!! Loop over all the particles
        Ptot = PREF*data(1,IT)         ! the total momentum
        Etot = SQRT(Ptot*Ptot + AM*AM) ! the total energy
        phi   = data(5,IT)*0.001D0     ! angle phi
        theta = data(3,IT)*0.001D0     ! angle theta

C---------Convert the coordinate system from Zgobi to DTA
        z_dta = data(6,IT)*cos(phi)*cos(theta)
        x_dta = data(2,IT)
        y_dta = data(4,IT)

        px_dta = Ptot*cos(phi)*sin(theta)
        py_dta = Ptot*sin(phi)

        pt_dta = Etot
        t_dta = data(7,IT)        ! time = F(7,I), but TAR = F(7,I)*1.D5
        s_dta = data(6,IT)        ! displacement = F(6,I)

C---------Scale the coordinated from DTA to dimensionless
        X_scl = x_dta/scl_l0
        Y_scl = y_dta/scl_l0
        Z_scl = z_dta/scl_l0

        PX_scl = px_dta/PREF
        PY_scl = py_dta/PREF
        PT_scl = pt_dta/PREF

        T_scl = -t_dta*scl_w0
        S_scl =  s_dta/scl_l0

C---------Start the symplectic dfift-kick-drift integrator
        pt_gb2 = PT_scl*PT_scl - gb0*gb0
        NUM_STEP = 1

        DO 999 WHILE (NUM_STEP .LE. MAX_STEP)
C----------DRIFT half stepsize for the first drift
          IF (NUM_STEP .EQ. 1) THEN
            ps = SQRT(pt_gb2 - PX_scl*PX_scl - PY_scl*PY_scl)
            portion = hfPAS_scl/ps

            X_scl = X_scl + portion*PX_scl
            Y_scl = Y_scl + portion*PY_scl
            DTscl = portion*PT_scl
            T_scl = T_scl - DTscl

            deltaZ = hfPAS_scl*
     >             SQRT(1-(PX_scl*PX_scl+PY_scl*PY_scl)/(ps*ps))
            Z_scl = Z_scl + deltaZ
            S_scl = S_scl + hfPAS_scl
          ENDIF

C----------KICK----------------------------------------------
          PX_scl = PX_scl - eb2*X_scl
          PY_scl = PY_scl + eb2*Y_scl

C----------DRIFT the whole step size except for the last step
          ps = SQRT(pt_gb2 - PX_scl*PX_scl - PY_scl*PY_scl)
          IF (NUM_STEP .EQ. MAX_STEP) THEN
            portion = hfPAS_scl/ps
            deltaZ = hfPAS_scl*
     >             SQRT(1-(PX_scl*PX_scl+PY_scl*PY_scl)/(ps*ps))
            S_scl = S_scl + hfPAS_scl
          ELSE
            portion = PAS_scl/ps
            deltaZ = PAS_scl*
     >             SQRT(1-(PX_scl*PX_scl+PY_scl*PY_scl)/(ps*ps))
            S_scl = S_scl + PAS_scl
          ENDIF

          X_scl = X_scl + portion*PX_scl
          Y_scl = Y_scl + portion*PY_scl
          DTscl = portion*PT_scl
          T_scl = T_scl - DTscl
          Z_scl = Z_scl + deltaZ

          NUM_STEP = NUM_STEP + 1
 999    CONTINUE

C----------Revert scaling coordinates from dimensionless to DTA
        x_dta = X_scl*scl_l0
        y_dta = Y_scl*scl_l0
        z_dta = Z_scl*scl_l0

        px_dta = PX_scl*PREF
        py_dta = PY_scl*PREF

        t_dta  = -T_scl/scl_w0
        s_dta  = S_scl*scl_l0

C----------Revert coordinates from DTA to Zgoubi
C----------Assuming phi and theta are both quite small
        Pxz = SQRT(Ptot*Ptot - py_dta*py_dta) ! the projected momentum on the xz-plane
        theta = ASIN(px_dta/Pxz)  ! theta = arcsin(px/p_y)
        phi = ASIN(py_dta/Ptot)   ! phi = arcsin(py/P)

        X_axis(IT) = z_dta
        data(2,IT) = x_dta
        data(3,IT) = theta*1000.D0
        data(4,IT) = y_dta
        data(5,IT) = phi*1000.D0
        data(6,IT) = s_dta
        data(7,IT) = t_dta

      ENDDO

      RETURN
      END
