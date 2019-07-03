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
      SUBROUTINE KDK_INTEGR(MAX_STEP, PAS_scl, gb0, eb2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C     --------------------------------------------------------------------
C     Kick(L/2)Drift(L)Kick(L/2) motion integrator for magnetic quadrupole
C     ALL the variables are scaled
C     --------------------------------------------------------------------

      INCLUDE "C.TRAJ_SCL.H" ! COMMON/TRAJ_SCL/ X_scl, PX_scl, Y_scl, PY_scl, T_scl, PT_scl, S_scl, Z_scl

      pt_gb2 = PT_scl*PT_scl - gb0*gb0
      NUM_STEP = 0

      DO 999 WHILE (NUM_STEP .LT. MAX_STEP)
C----------KICK
        PX_scl = PX_scl - eb2*X_scl
        PY_scl = PY_scl + eb2*Y_scl
C----------DRIFT
        ps = SQRT(pt_gb2 - PX_scl*PX_scl - PY_scl*PY_scl)
        portion = PAS_scl/ps
        X_scl = X_scl + portion*PX_scl
        Y_scl = Y_scl + portion*PY_scl

        DTscl = portion*PT_scl
        T_scl = T_scl - DTscl

        deltaZ = PAS_scl*SQRT(1-(PX_scl*PX_scl+PY_scl*PY_scl)/(ps*ps))
        Z_scl = Z_scl + deltaZ
        S_scl = S_scl + PAS_scl
C----------KICK
        PX_scl = PX_scl - eb2*X_scl
        PY_scl = PY_scl + eb2*Y_scl

        NUM_STEP = NUM_STEP + 1
 999  CONTINUE

c      WRITE(*, 101) NUM_STEP, Z_scl*100.D0, MAX_STEP
c 101  FORMAT(/, 6X, 'After ', I8, ' steps, X = ', E16.8, ' cm')

      RETURN
      END
