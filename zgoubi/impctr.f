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
      SUBROUTINE IMPCTR(IUNIT,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXV=40) 
      COMMON/CONTR/VAT(MXV),XI(MXV)
      COMMON /VAR/ X(3*MXV),P(MXV)
      COMMON/VARY/NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,7)
 
        WRITE(IUNIT,500)
500     FORMAT(' STATUS OF CONSTRAINTS')
        WRITE(IUNIT,600)
600     FORMAT(
     >  ' TYPE  I   J  LMNT#      DESIRED         WEIGHT        ',
     >  '  REACHED         KI2         *  Parameter(s) ')
        DO 2 I=1,NC
          XI2=((VAT(I)-V(I))/W(I))**2/F
          NPRM1 = NINT(CPAR(I,1)) + 1
          WRITE(IUNIT,700) IC(I),I1(I),I2(I),I3(I),V(I),W(I),VAT(I),XI2, 
     >    NINT(CPAR(I,1)),(CPAR(I,JJ),JJ=2,NPRM1)
700       FORMAT(3I4,I6,5X,1P,G12.5,4X,G11.4,3X,G14.7,2X,G11.4,2X,
     >    '   * ',I2,' : ',6(G9.1,'/'))
2      CONTINUE
 
      RETURN
      END
