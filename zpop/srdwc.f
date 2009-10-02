C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  François Méot
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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE SRDWC(MNU,HZ,R,
     >                               DW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MSAM=1000)
      DIMENSION DW(MSAM,*)

      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH
      INCLUDE 'MXSTEP.H'
      COMMON/SRTAB/ EFD(MXSTEP,3),D3W(MSAM,9),NPTS

      DIMENSION SR(3), SI(3)

      DIMENSION ANU(MSAM), DNU(MSAM)
      COMMON/FREQ/ ANU, DNU

      CMU = 4.D-7 * PI * CL
      R2CMU =   2.D0 * R*R / CMU

      DO 1 IN = 1, MNU
        OMG = DPI * ANU(IN) * HZ                       ! rd/s
        DT = .5D0 * EFD(1,1)                       ! dt (s)
C        DT = EFD(1,1)                       ! dt (s)
        OMT = 0.D0
C------- Y component
        SR(2) = EFD(1,2) * DT                      ! Ey.dt (V.s)
        SI(2) = 0.D0  
C------- Z component
        SR(3) = EFD(1,3) * DT                      ! Ez.dt
        SI(3) = 0.D0

        DO 2 I = 2, NPTS - 1
          DT = 0.5D0 * (EFD(I-1,1) + EFD(I,1))           ! dt (s)
          OMT = OMT + OMG * DT                           ! omga.t (rad)
          IF( OMT .GE. DPI ) OMT = OMT -DPI 
          COMT = COS(OMT)
          SIMT = SIN(OMT)

          FC = EFD(I,2) * DT

C--------- Y component
          SR(2) = SR(2) + FC * COMT
          SI(2) = SI(2) - FC * SIMT
C--------- Z component
          FC = EFD(I,3) * DT
          SR(3) = SR(3) + FC * COMT
          SI(3) = SI(3) - FC * SIMT
 2      CONTINUE

        DT = 0.5D0 * EFD(NPTS - 1,1)
C        DT = EFD(NPTS - 1,1)
        OMT = OMT + OMG * DT 
        IF( OMT .GE. DPI ) OMT = OMT -DPI 
        COMT = COS(OMT)
        SIMT = SIN(OMT)

        FC = EFD(NPTS,2) * DT

C------- Y component
        SR(2) = SR(2) + FC * COMT
        SI(2) = SI(2) - FC * SIMT
C------- Z component
        FC = EFD(NPTS,3) * DT
        SR(3) = SR(3) + FC * COMT
        SI(3) = SI(3) - FC * SIMT

        DW(IN,1) = ANU(IN)

C------- Energy spectrum dW/dNu.dOmga (J/Hz.srd)
C        in the y  and z field  components.
        DW(IN,2) = R2CMU *  ( SR(2)*SR(2) + SI(2)*SI(2) )
        DW(IN,3) = R2CMU *  ( SR(3)*SR(3) + SI(3)*SI(3) )

 1    CONTINUE

      RETURN
      END
