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
      SUBROUTINE MKM_INTEGR(MAX_STEP, PAS, b0g0, eb2, AM, PREF, IMAX,
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

      DOUBLE PRECISION theta, phi, x, y, z, s, t, px, py, pt, S_scl
      DOUBLE PRECISION Z_scl, X_scl, PX_scl, Y_scl, PY_scl, T_scl
      DOUBLE PRECISION PT_scl, X0, PX0, Y0, PY0, X1, PX1, Y1, PY1
      DOUBLE PRECISION delX, delY, delZ, hfPAS, PAS2, delP, delT, Ptot
      DOUBLE PRECISION Etot, P1, P2, P3, P1m1, tP3, Ps, Psm1, Pxy2, Pxz
      DOUBLE PRECISION kp, kl, ckl, skl, chkl, shkl
      DOUBLE PRECISION tkl, c2kl, s2kl, ch2kl, sh2kl
      DOUBLE PRECISION const1, const2, const3, kappa
      DOUBLE PRECISION MQ11, MQ12, MQ21, MQ22, MQ31, MQ32, MQ41, MQ42
      DOUBLE PRECISION MQ51, MQ52, MQ53, MQ54, MQ55, MQ56
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
      hfPAS = 0.5D0*PAS                ! 0.5L/L0
      PAS2  = PAS*PAS                  ! PAS = L/L0, scaled

      DO CONCURRENT (IT = 1:IMAX)  !!! Loop over all the particles
        Ptot = PREF*data(1,IT)         ! the total momentum
        Etot = SQRT(Ptot*Ptot + AM*AM) ! the total energy
        phi   = data(5,IT)*0.001D0     ! angle phi
        theta = data(3,IT)*0.001D0     ! angle theta

C---------Convert the coordinates from Zgobi notation to DTA notation---
        z  = data(6,IT)*cos(phi)*cos(theta)
        x  = data(2,IT)
        y  = data(4,IT)
        px = Ptot*cos(phi)*sin(theta)
        py = Ptot*sin(phi)
        pt = Etot
        t  = data(7,IT)        ! time = F(7,I), but TAR = F(7,I)*1.D5
        s  = data(6,IT)        ! displacement = F(6,I)

C---------Scale the coordinated from DTA to dimensionless---------------
        X_scl = x/l0
        Y_scl = y/l0
        Z_scl = z/l0
        PX_scl = px/PREF
        PY_scl = py/PREF
        PT_scl = pt/PREF
        T_scl = -t*w0
        S_scl =  s/l0

C---------Compute constants for the symplectic matrix-kick-matrix integrator
C---------These constants do not change for a particle during MKM-------
        P2 = PT_scl*PT_scl - b0g0*b0g0      ! P^2 = P_T^2 - (1/(beta0*gamma0))^2
        P1 = SQRT(P2)                       ! P
        P3 = P1**3                          ! P^3
        P1m1 = 1.D0/P1                      ! P^(-1)
        tP3  = 2.D0*P3                      ! 2P^3
        kappa  = SQRT(eb2/P1)               ! K
        const1 = 0.125D0*PT_scl/P3/kappa    ! P_T/2P/P^2/4K = 0.125*P_T/P^3/K
        const2 = 0.125D0*PT_scl*kappa/P1    ! P_T*K^2/2P/4K = 0.125*P_T*K/P
        const3 = 0.250D0*PT_scl/P2          ! P_T*K/2P/P/2K = 0.250*P_T/P^2

        kp = kappa*P1               ! KP is a constant for a particle
        kl = kappa*hfPAS            ! KL is a constant for a particle
                                    ! Note that M_Q(L/2), half of L=PAS
        ckl = cos(kl)
        skl = sin(kl)
        chkl = cosh(kl)
        shkl = sinh(kl)

        tkl = kappa*PAS               ! 2KL, where L/2 is for M_Q operation
        s2kl  = sin(tkl)              ! sin(2KL)
        c2kl  = cos(tkl)              ! cos(2KL)
        sh2kl = sinh(tkl)             ! sinh(2KL)
        ch2kl = cosh(tkl)             ! cosh(2KL)

        MQ11 = ckl                    ! MQ?? are transition maxtrix elements
        MQ12 = skl/kp
        MQ21 = ckl
        MQ22 = -kp*skl
        MQ31 = chkl
        MQ32 = shkl/kp
        MQ41 = chkl
        MQ42 = shkl*kp

        MQ51 = const1*(tkl+s2kl)   ! (P_T/2P^3)*[2KL+sin(2KL)]/4K
        MQ52 = const1*(sh2kl+tkl)  ! (P_T/2P^3)*[sinh(2KL)+2KL]/4K
        MQ53 = const2*(tkl-s2kl)   ! (P_T*K^2/2P)*[2KL-sin(2KL)]/4K
        MQ54 = const2*(sh2kl-tkl)  ! (P_T*K^2/2P)[sinh(2KL)-2KL]/4K
        MQ55 = const3*(c2kl-1.D0)  ! (P_T*K/2P^2)*[cos(2KL)-1]/2K
        MQ56 = const3*(ch2kl-1.D0) ! (P_T*K/2P^2)*[cosh(2KL)-1]/2K

C---------Start the symplectic matrix-kick-matrix integrator------------
        NUM_STEP = 1
        DO 999 WHILE (NUM_STEP .LE. MAX_STEP)
C----------M_Q(L/2)-----------------------------------------------------
          X0  = X_scl
          Y0  = Y_scl
          PX0 = PX_scl
          PY0 = PY_scl

          X_scl  =  X0*MQ11 + PX0*MQ12
          PX_scl = PX0*MQ21 +  X0*MQ22
          Y_scl  =  Y0*MQ31 + PY0*MQ32
          PY_scl = PY0*MQ41 +  Y0*MQ42
          T_scl  = T_scl -
     >            (PX0*PX0*MQ51 + PY0*PY0*MQ52
     >            + X0* X0*MQ53 +  Y0* Y0*MQ54
     >            + X0*PX0*MQ55 +  Y0*PY0*MQ56)

C----------K_Q(L)-------------------------------------------------------
          Pxy2 = PX_scl*PX_scl + PY_scl*PY_scl
          Ps   = SQRT(P2 - Pxy2)
          Psm1 = 1.D0/Ps
          delP = (Psm1 - P1m1)*PAS
          delT = (Psm1 - Pxy2/tP3)*PAS

          X_scl = X_scl + delP*PX_scl
          Y_scl = Y_scl + delP*PY_scl
          T_scl = T_scl - delT*PT_scl

C----------M_Q(L/2)-----------------------------------------------------
          X1  = X_scl
          Y1  = Y_scl
          PX1 = PX_scl
          PY1 = PY_scl

          X_scl  =  X1*MQ11 + PX1*MQ12
          PX_scl = PX1*MQ21 +  X1*MQ22
          Y_scl  =  Y1*MQ31 + PY1*MQ32
          PY_scl = PY1*MQ41 +  Y1*MQ42
          T_scl  = T_scl -
     >            (PX1*PX1*MQ51 + PY1*PY1*MQ52
     >            + X1* X1*MQ53 +  Y1* Y1*MQ54
     >            + X1*PX1*MQ55 +  Y1*PY1*MQ56)

C-----------Update Z_scl and S_scl after one cycle of MKM---------------
          delX = X_scl - X0
          delY = Y_scl - Y0
          delZ = SQRT(PAS2 - delX**2 - delY**2)
          Z_scl = Z_scl + delZ
          S_scl = S_scl + PAS

          NUM_STEP = NUM_STEP + 1
 999    CONTINUE

C----------Revert scaling coordinates from dimensionless to DTA---------
        x  = X_scl*l0
        y  = Y_scl*l0
        z  = Z_scl*l0
        px = PX_scl*PREF
        py = PY_scl*PREF
        t  = -T_scl/w0
        s  = S_scl*l0

C----------Revert coordinates from DTA notation to Zgoubi---------------
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
