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
      SUBROUTINE ZEROF(IT,XL,BTA,ALP,OME,C1,R,EPS,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ** VALEUR DE T QUI EST LE ZERO DE Y(T)=XL
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M

      T = XL/(BTA*CL)
      Y = -R*COS(OME*T+EPS) - (ALP*T  -C1)/OME
      IF( (Y-XL)*(Y-XL) .LT. 1.D-6 ) RETURN
      N=0
 1    CONTINUE
      N = N+1
      T = ( C1 - ( R*COS(OME*T+EPS) + XL ) *OME ) / ALP
      Y = -R*COS(OME*T+EPS) - (ALP*T  -C1)/OME
      IF( (Y-XL)*(Y-XL) .LT. 1.D-6 ) RETURN
      IF(N .GT. 200) GOTO 98
      GOTO 1

 98   IF(NRES.GT.0) WRITE(NRES,*) ' TRAJECTOIRE ',IT
     > ,', Convergence problem in sbr ZEROF :',N,'  iterations'
      RETURN
      END
