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
      SUBROUTINE CONVERT_COORD(PREF, Ptot, Etot)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
C     --------------------------------------------------------------
C     Convert the coordinates of a particle from the Zgoubi notation 
C     to the notation by DTA (Dan T. Abell)
C     --------------------------------------------------------------

      INCLUDE "C.TRAJ.H"     ! COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      INCLUDE "C.TRAJ_DTA.H" ! COMMON/TRAJ_DTA/ x_dta, px_dta, y_dta, py_dta, t_dta, pt_dta, s_dta, z_dta
      INCLUDE "C.TRAJ_SCL.H" ! COMMON/TRAJ_SCL/ X_scl, PX_scl, Y_scl, PY_scl, T_scl, PT_scl, S_scl, Z_scl
      INCLUDE "SCALE_KDK.H"  ! PARAMETER (scl_l0, scl_w0)
      
C---------Convert the coordinate system 
     
      z_dta = X
      x_dta = Y
      y_dta = Z
      
      px_dta = Ptot*cos(P)*sin(T)  ! p_x = Pcos(phi)sin(theta)      
      py_dta = Ptot*sin(P)         ! p_y = Psin(phi)
      
      pt_dta = Etot
      t_dta = TAR/1.D5             ! time = F(7,I), but TAR = F(7,I)*1.D5
      s_dta = SAR                  ! displacement = F(6,I)
      
c      WRITE(*, 100) 'Before kdk, converted to DTA''s notation:'
c      WRITE(*, 101) x_dta, px_dta, y_dta, py_dta, 
c     >  t_dta, pt_dta, s_dta, z_dta
      
c 100  FORMAT(/, 10X, A)
c 101  FORMAT(15X, 'x = ', E16.8, ' cm, px = ', E16.8, ' MeV/c',
c     >/, 15X, 'y = ', E16.8, ' cm, py = ', E16.8, ' MeV/c',
c     >/, 15X, 't = ', E16.8, ' s , pt = ', E16.8, ' MeV/c2',
c     >/, 15X, 's = ', E16.8, ' cm, z  = ', E16.8, ' cm')
     
C---------Scale the coordinated from DTA to dimensionless
      
      X_scl = x_dta/scl_l0
      Y_scl = y_dta/scl_l0
      Z_scl = z_dta/scl_l0
      
      PX_scl = px_dta/PREF
      PY_scl = py_dta/PREF
      PT_scl = pt_dta/PREF
      
      T_scl = -t_dta*scl_w0
      S_scl = s_dta/scl_l0
      
c      WRITE(*, 100) 'Before kdk, the scaled variables:'
c      WRITE(*, 102) X_scl, PX_scl, Y_scl, PY_scl,
c     >  T_scl, PT_scl, S_scl, Z_scl
c 102  FORMAT(15X, 'X = ', E16.8, ', PX = ', E16.8,
c     >/, 15X, 'Y = ', E16.8, ', PY = ', E16.8,
c     >/, 15X, 'T = ', E16.8, ', PT = ', E16.8,
c     >/, 15X, 'S = ', E16.8, ',  Z = ', E16.8)
 
      RETURN
      END
      
      