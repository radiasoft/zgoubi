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
      SUBROUTINE FBORNE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER * 81 TAMP
      COMMON/EDIT/TAMP
      PARAMETER (MXV=40) 
      COMMON/MIMA/ DX(MXV),XMI(MXV),XMA(MXV)
      COMMON /VAR/ X(3*MXV),P(MXV)
      COMMON/VARY/NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,7)
      DO I=1,NV
C FM, Dec.2002, for fit of longitudinal ellipses out of pi-collect channel
C            P(I)=ABS(DX(I)*X(I)/10.D0) +.01  replaced by
C            P(I)=ABS(DX(I)*X(I)/10.D0) 
C---------------------------
C           P(I)=ABS(DX(I)*X(I)/10.D0) +.01D0  
           P(I)=(XMA(I)-XMI(I))/100.D0  
           K=I+NV
           J=K+NV
           KL=XCOU(I)
           IF(KL .EQ. 0) THEN
             X(K)=XMI(I)
             X(J)=XMA(I)
           ELSE
             IF(KL.LT.0) THEN
               X(K)=-XMA(I)
               X(J)=-XMI(I)
             ELSE
               X(K)=XMI(I)
               X(J)=XMA(I)
             ENDIF
           ENDIF
      ENDDO
      RETURN
      END
