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
C  USA
C  -------
      FUNCTION AMINI(FNCT,CX,CY,V,VMI,VMA,NC,NV,RNU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(MXV=10)
      PARAMETER(MXC=400)
      DIMENSION V(MXV),CX(MXC),CY(MXC),FCT(MXC)
      DIMENSION RNU(MXV),VMI(MXV),VMA(MXV)
      DIMENSION V1(MXV),V2(MXV),AL1(MXV),AL2(MXV)

      EXTERNAL FNCT

      DATA TO/ .6180339885D0/

      DO 10 I=1,NV
        AA=(VMI(I)-V(I))/RNU(I)
        BB=(VMA(I)-V(I))/RNU(I)
        AL1(I)=DMIN1(AA,BB)
10      AL2(I)=DMAX1(AA,BB)
      A0=DMAX1(AL1(1),AL1(2),AL1(3))
      B0=DMIN1(AL2(1),AL2(2),AL2(3))
      DU=(1.D0-TO)*(B0-A0)
      U1=A0+DU
      U2=B0-DU
      DO 15 I=1,NV
        V1(I)=V(I)+U1*RNU(I)
15      V2(I)=V(I)+U2*RNU(I) 
      CALL FUNK(FNCT,V1,CX(MXC-1),CX(MXC),NC,NV,CX,FCT)
      SQ1=SQUARE(CY,NC,FCT)
      CALL FUNK(FNCT,V2,CX(MXC-1),CX(MXC),NC,NV,CX,FCT)
      SQ2=SQUARE(CY,NC,FCT)
      AMINI=DOREE(FNCT,A0,B0,U1,U2,V,CX,CY,NC,NV,RNU,SQ1,SQ2)
      RETURN
      END
