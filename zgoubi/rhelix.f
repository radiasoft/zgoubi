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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE RHELIX(ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL 

      READ(NDAT,*) IA
      A(NOEL,1) = IA
C----- Length (cm), Twist-pitch=distance for 2pi B rotation (cm), BO (kG), initial angle of \vecB wrt vertical axis (rad)
      READ(NDAT,*) (A(NOEL,I),I=10,13)
C----- Kird, resol. Anal :  kird=0,  resol= 2(default) or 4.   Num :  kird=3, resol = xpas/meshSize    
      READ(NDAT,*) A(NOEL,20), A(NOEL,21)

      ND = 30
      CALL STPSIZ(NDAT,NOEL,ND,
     >                         A)

C----- KP,XCE, YCE, ALE
      ND1=ND+10
      READ(NDAT,*) (A(NOEL,I),I=ND1,ND1+3)
      RETURN
      END
