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
C  François Meot <fmeot@bnl.gov>
C  Brookhaven National Laboratory
C  C-AD, Bldg 911
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE FBORNE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      PARAMETER (MXV=60)
      INCLUDE "C.MIMA.H"     ! COMMON/MIMA/ DX(MXV),XMI(MXV),XMA(MXV)
      INCLUDE "C.VAR.H"     ! COMMON /VAR/ X(3*MXV),P(MXV)
      INCLUDE "C.VARY.H"  ! COMMON/VARY/ NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
                          !     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,27)
      DO I=1,NV
C FM, Dec.2002, for fit of longitudinal ellipses out of pi-collect channel
C            P(I)=ABS(DX(I)*X(I)/10.D0) +.01  replaced by
C            P(I)=ABS(DX(I)*X(I)/10.D0)
C---------------------------
C           P(I)=ABS(DX(I)*X(I)/10.D0) +.01D0
           P(I)=(XMA(I)-XMI(I))/100.D0
           K=I+NV
           J=K+NV
           KL=INT(XCOU(I))
C FM Jan 2015. This was a mistake. The coupled variable is just forced to
C the value reached by the varied variable (in sbr rempli)
c           IF(KL .EQ. 0) THEN
             X(K)=XMI(I)
             X(J)=XMA(I)
c           ELSE
c             IF(KL.LT.0) THEN
c               X(K)=-XMA(I)
c               X(J)=-XMI(I)
c             ELSE
c               X(K)=XMI(I)
c               X(J)=XMA(I)
c             ENDIF
c           ENDIF
      ENDDO
      RETURN
      END
