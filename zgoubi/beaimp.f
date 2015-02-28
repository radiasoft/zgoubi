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
C  Brookhaven National Laboratory    
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE BEAIMP(F0,PHY,PHZ) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/BEAM/ FI(6,6)
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      DIMENSION F0(6,*)

      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      IF(NRES .LT. 0) RETURN
      WRITE(NRES,103) 'INITIAL'
 103  FORMAT(//,18X,'BEAM  MATRIX (beta/alpha/alpha/gamma, D,D''), 
     >       ',A,/)
      WRITE(NRES,104) (( FI(IA,IB) , IB=1,6) , IA=1,6)
 104  FORMAT(6X,1P,6G16.6)
      WRITE(NRES,103) 'FINAL' 
      WRITE(NRES,104) (( F0(IA,IB) , IB=1,6) , IA=1,6)
      WRITE(NRES,FMT='(/,18X,''Betatron phase advances (fractional),  ''
     >,''phi_y/2pi,'',''  phi_z/2pi :''
     >,//,18X,1P,2(5X,E14.6))') PHY/(2.D0*PI), PHZ/(2.D0*PI)
      RETURN
      END
