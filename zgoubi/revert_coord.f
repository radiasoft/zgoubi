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
C  Upton, NY, 11973
C  -------
      SUBROUTINE REVERT_COORD(PREF, Ptot)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
C     --------------------------------------------------------------
C     Revert the coordinates of a particle from
C     DTA's notation to the Zgoubi notation
C     --------------------------------------------------------------

      INCLUDE "C.TRAJ.H"     ! COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      INCLUDE "C.TRAJ_DTA.H" ! COMMON/TRAJ_DTA/ x_dta, px_dta, y_dta, py_dta, t_dta, pt_dta, s_dta, z_dta
      INCLUDE "C.TRAJ_SCL.H" ! COMMON/TRAJ_SCL/ X_scl, PX_scl, Y_scl, PY_scl, T_scl, PT_scl, S_scl, Z_scl
      INCLUDE "SCALE_KDK.H"  ! PARAMETER (scl_l0, scl_w0)

C----------revert scaling coordinates from dimensionless to DTA
      scl_p0 = PREF

      x_dta = X_scl*scl_l0
      y_dta = Y_scl*scl_l0
      z_dta = Z_scl*scl_l0

      px_dta = PX_scl*PREF
      py_dta = PY_scl*PREF

      t_dta  = -T_scl/scl_w0
      s_dta  = S_scl*scl_l0
      
C----------Revert coordinates from DTA to Zgoubi
C----------Assuming phi and theta are both small
      Y = x_dta
      Z = y_dta
      X = z_dta

      Pxz = SQRT(Ptot*Ptot - py_dta*py_dta) ! the projected momentum on the xz-plane
      T = ASIN(px_dta/Pxz)  ! theta = arcsin(px/p_y)
      P = ASIN(py_dta/Ptot) ! phi = arcsin(py/P)

      SAR = s_dta        ! displacement
      TAR = t_dta*1.D5   ! time

c      WRITE(*, 100) 'After kdk integration in Dan''s notation:'
c      WRITE(*, 101) x_dta, px_dta, y_dta, py_dta,
c     >  t_dta, pt_dta, s_dta, z_dta

c 100  FORMAT(/, 10X, A)
c 101  FORMAT(15X, 'x = ', E16.8, ' cm, px = ', E16.8, ' MeV/c',
c     >/, 15X, 'y = ', E16.8, ' cm, py = ', E16.8, ' MeV/c',
c     >/, 15X, 't = ', E16.8, ' s , pt = ', E16.8, ' MeV/c2',
c     >/, 15X, 's = ', E16.8, ' cm, z  = ', E16.8, ' cm')

c      WRITE(*, 102) scl_w0, PREF, Ptot, Pxz
c 102  FORMAT(/, 10X, 'w0 = ', E16.8, 4X, 'PREF = ', E16.8,
c     >  /, 10X, 'Ptot = ', E16.8, 4X, 'Pxz = ', E16.8)
      
      RETURN
      END
      
      
