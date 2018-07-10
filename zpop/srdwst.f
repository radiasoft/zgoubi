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
C  François Méot <fmeot@bnl.gov>
C  Brookhaven National Laboratory                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE SRDWST(MNU,JY,KZ,DW,HZ,DOM,DPH,DPS,
     >                                              D3W,SDW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MSAM=1000)
      DIMENSION DW(MSAM,*), SDW(*), D3W(MSAM,*)
C     ------------------------------------------------
C     Sum and store spectrum, fragment after fragment.
C     ------------------------------------------------
      COMMON/CDF/ IES,IORDRE,LCHA,LIST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN

      DIMENSION ANU(MSAM), DNU(MSAM)
      COMMON/FREQ/ ANU, DNU

C------ DW(IN,2-3) = Energy spectrum dW_y,z/dNu.dOmga (J/Hz.srd) at current Phi, Psi
      IN = 1
      DNUI = DNU(IN)*0.5D0
      Y = DW(IN,2) *DOM
      Z = DW(IN,3) *DOM
      SDW(2) = Y * DNUI                         
      SDW(3) = Z * DNUI
      D3W(IN,2) = D3W(IN,2) + Y                    
      D3W(IN,3) = D3W(IN,3) + Z
      DO 1 IN=2,MNU-1
        DNUI = DNU(IN)            !erreur corrction may 99 *0.5D0
        Y = DW(IN,2) *DOM 
        Z = DW(IN,3) *DOM 
        SDW(2) = SDW(2) + Y * DNUI                         ! sum over Nu
        SDW(3) = SDW(3) + Z * DNUI
        D3W(IN,2) = D3W(IN,2) + Y         ! dWy,z over Y strip
        D3W(IN,3) = D3W(IN,3) + Z
 1    CONTINUE
      IN = MNU
      DNUI = DNU(IN)*0.5D0
      Y = DW(IN,2) *DOM 
      Z = DW(IN,3) *DOM 
      SDW(2) = SDW(2) + Y * DNUI
      SDW(3) = SDW(3) + Z * DNUI
C----- End up with Spectral density at current fragmnt
      D3W(IN,2) = D3W(IN,2) + Y           ! dW/dNu (J/Hz)
      D3W(IN,3) = D3W(IN,3) + Z

C----- Total energy (J) in the current fragment :
      SDW(2) = SDW(2) * HZ
      SDW(3) = SDW(3) * HZ

C----- Angular density at current Phi (Y) angle [integrated over Psi (Z)]
      D3W(JY,5) = D3W(JY,5) + SDW(2) / DPH               ! dW/dPhi (J/rad)
      D3W(JY,6) = D3W(JY,6) + SDW(3) / DPH
C----- Angular density at current Psi (Z) angle [integrated over Phi (Y)]
      D3W(KZ,8) = D3W(KZ,8) + SDW(2) / DPS               ! dW/dPsi (J/rad)
      D3W(KZ,9) = D3W(KZ,9) + SDW(3) / DPS

      RETURN
      END
