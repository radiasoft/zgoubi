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
C  Brookhaven National Laboratory                                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
C                        
      FUNCTION DOREE(FNCT,A0,B0,U1,U2,V,CX,CY,NC,NV,RNU,G1,G2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(MXV=10)
      PARAMETER(MXC=400)
      DIMENSION V(MXV),CX(MXC),CY(MXC),FCT(MXC)
      DIMENSION RNU(MXV)
      DIMENSION DV(MXV)
      EXTERNAL FNCT

      DATA TO/.6180339885D0/

60    CONTINUE
        IF    (G1.LT. G2) THEN
          B0=U2
          IF(ABS((A0-B0)/B0).GT. 1.D-5) THEN
            U2=U1
            G2=G1
            U1=A0+(1.D0-TO)*(B0-A0)
            DO 10 I=1,NV
10            DV(I)=V(I)+U1*RNU(I)
            CALL FUNK(FNCT,DV,CX(MXC-1),CX(MXC),NC,NV,CX,FCT)
            G1=SQUARE(CY,NC,FCT)
            GOTO 60
          ENDIF
        ELSEIF(G1.GE. G2) THEN
          A0=U1
          IF(ABS((A0-B0)/B0).GT. 1.D-5) THEN
            U1=U2
            G1=G2
            U2=B0-(1.D0-TO)*(B0-A0)
            DO 20 I=1,NV
20            DV(I)=V(I)+U2*RNU(I)
            CALL FUNK(FNCT,DV,CX(MXC-1),CX(MXC),NC,NV,CX,FCT)
            G2=SQUARE(CY,NC,FCT)
            GOTO 60
          ENDIF
        ENDIF
C
      DOREE=(A0+B0)*.5D0
      RETURN
      END
