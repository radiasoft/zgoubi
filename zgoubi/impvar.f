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
      SUBROUTINE IMPVAR(IUNIT,NI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXV=40) 
      COMMON/CONTR/VAT(MXV),XI(MXV)
      COMMON /VAR/ X(3*MXV),P(MXV)
      COMMON/VARY/NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
     >IC(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,7)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
 
        IF(IUNIT .EQ. 7) WRITE(IUNIT,100)
100     FORMAT('1')
        WRITE(IUNIT,200) NI
200     FORMAT(/,' STATUS OF VARIABLES  (Iteration #',I6,')')
        WRITE(IUNIT,300)
300     FORMAT(
     >  ' LMNT  VAR  PARAM  MINIMUM     INITIAL        FINAL        ',
C------   IR(I)  I   IS(I)   X(K)        XI(I)         X(I)  
     >  'MAXIMUM     STEP' )
C----     X(J)        P(I)
 
      DO 1 I=1,NV
        K=I+NV
        J=K+NV
        WRITE(IUNIT,400) IR(I),I,IS(I),X(K),XI(I),
     >  A(IR(I),IS(I)),X(J),P(I)
 400    FORMAT(1P, 
     >  2X,I3,3X,I2,4X,I3,2(2X,G10.3),2X,G15.8,2(1X,G10.3))
        IF(XCOU(I).NE.0.D0) THEN
          KL=XCOU(I)
          KP=NINT((100.D0*XCOU(I)-100.D0*KL))
          WRITE(IUNIT,400) KL,I,KP,
     >            A(IR(I),IS(I)),XI(I),X(I),X(J),P(I)
        ENDIF
 1    CONTINUE
      RETURN
      END
