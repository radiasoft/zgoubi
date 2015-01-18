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
      SUBROUTINE FITARR(IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXV=60)
      COMMON/VARY/NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,27)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      IER=0
      IF(NC.LT.2) RETURN
      DO 4 J=2,NC
        II1=J-1
        DO 5 I=J,NC
          IF(I3(I).GE.I3(II1))   GOTO 5
          IC0=IC(I)
          I10=I1(I)
          I20=I2(I)
          I30=I3(I)
          V0=V(I)
          W0=W(I)
          IC(I)=IC(II1)
          I1(I)=I1(II1)
          I2(I)=I2(II1)
          I3(I)=I3(II1)
          V(I)=V(II1)
          W(I)=W(II1)
          IC(II1)=IC0
          I2(II1)=I20
          I1(II1)=I10
          I3(II1)=I30
          V(II1)=V0
          W(II1)=W0
 5      CONTINUE
 4    CONTINUE
      IF(I3(NC).LE.NB)   GOTO 6
      IER=1
      WRITE(6,100)
 100  FORMAT(' CONTRAINTE HORS LIMITES=DEROUTE')
 6    CONTINUE
      RETURN
      END
