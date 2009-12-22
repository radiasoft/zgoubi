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
      SUBROUTINE FBORNE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER * 81 TAMP
      COMMON/EDIT/TAMP
      COMMON/MIMA/ DX(40)
      PARAMETER (MXV=40) 
      COMMON /VAR/ X(3*MXV),P(MXV)
      COMMON/VARY/NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,7)
 
      IBEA=0
      DO 1 I=1,NV
         IF(IR(I) .EQ. 0) THEN
            IBEA=IBEA+1
         ELSE
C FM, Dec.2002, for fit of longitudinal ellipses out of pi-collect channel
C            P(I)=ABS(DX(I)*X(I)/10.D0) +.01  replaced by
C            P(I)=ABS(DX(I)*X(I)/10.D0) 
C---------------------------
            P(I)=ABS(DX(I)*X(I)/10.D0) +.01D0  
            K=I+NV
            J=K+NV
            KL=XCOU(I)
            IF(KL .NE. 0) THEN
               KP=NINT((1D3*XCOU(I)-1D3*KL))
               IF(KL.LT.0) THEN
C-------------- MINIMA
                 X(K)=0.D0
C-------------- MAXIMA
                 X(J)=X(I)+A(-KL,-KP)
               ELSE
                 X(K)=X(I)-ABS(X(I))*DX(I)
                 X(J)=X(I)+ABS(X(I))*DX(I)
               ENDIF
            ELSE
C-------------- MINIMA
               X(K)=X(I)-ABS(X(I))*DX(I)
C-------------- MAXIMA
               X(J)=X(I)+ABS(X(I))*DX(I)
            ENDIF
         ENDIF
 1    CONTINUE
      RETURN
      END
