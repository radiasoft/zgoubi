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
C  Upton, NY, 11973, USA
C  ------- 
c     propag.f
c     written by Frédéric Desforges, frederic.desforges@grenoble-inp.org
c     The propagation of the generalized twiss parameters uses a different 
c     mathematical process from the one which was presented in my report.
c     it uses the one given by Yun Luo in the paper entitled "Linear coupling
c     parametrization in the action-angle frame" of the PHYSICAL REVIEW
c     SPECIAL TOPICS - ACCELERATORS AND BEAMS, published 7 December 2004

      SUBROUTINE PROPAG(R,P,
     >                      F0,C,RPRM)
      IMPLICIT DOUBLE PRECISION (A-H,M-Z)
      DIMENSION F0(6,6),R(4,4),P(4,4)
      DIMENSION G(4,4),R11(2,2),R22INV(2,2),C(2,2),WORK(2)
      INTEGER I,J,N_ET,NTRANS,NRES,NTWISS,IPIV,INFO,DEBSTR,FINSTR
      DIMENSION IPIV(4)
      CHARACTER(300) BUFFER,LABEL,KEYWOR
      LOGICAL STRCON 
      DATA NU1, NU2 / 0.D0, 0.D0 /

C R is the 4x4 transport matrix from origin. See Luo : G yields coupled functions.
      G = MATMUL(R,P)

C     Propagation of the coupling parameters and the coupling strenght

      R11(1,1) = P(1,1)
      R11(1,2) = P(1,2)
      R11(2,1) = P(2,1)
      R11(2,2) = P(2,2)
      R22INV(1,1) = P(3,3)
      R22INV(1,2) = P(3,4)
      R22INV(2,1) = P(4,3)
      R22INV(2,2) = P(4,4)
      CALL DGETRF(2,2,R22INV,2,IPIV,INFO)
      CALL DGETRI(2,R22INV,2,IPIV,WORK,2,INFO)

      C = MATMUL(MATMUL(R11,C),R22INV)

      RPRM = 0.5*(SQRT(G(1,1)*G(2,2)-G(1,2)*G(2,1))+SQRT(G(3,3)*G(4,4)
     >-G(3,4)*G(4,3)))

      IF(RPRM .LE. 1. .AND. RPRM .GE. SQRT(2.)/2.) THEN
            CMOINS = 2*SQRT(ABS(1-1/RPRM**2))/(1+ABS(1-1/RPRM**2))*A
     >BS(NU1-NU2)
      ELSE IF(RPRM .LT. SQRT(2.)/2.) THEN
            RPRM = RPRM+2*(SQRT(2.)/2.-RPRM)                          ! According to the chosen convention 
            CMOINS = 2*SQRT(ABS(1-1/RPRM**2))/(1+ABS(1-1/RPRM**2))*A  ! RPRM cannot be lower than sqrt(2)/2
     >BS(NU2-NU1)                                                     ! For more detailed arguments look my internal
      ELSE
            CMOINS = 0.
      ENDIF


C     Propagation of the generalized Twiss' parameters and coupling parameters
      
      F0(1,1)  = (G(1,1)**2+G(1,2)**2) / (G(1,1)*G(2,2)-G(1,2)*G(2,1))
      F0(3,3)  = (G(3,3)**2+G(3,4)**2) / (G(3,3)*G(4,4)-G(3,4)*G(4,3))
      F0(1,2) = +( -(G(1,1)*G(2,1)+G(1,2)*G(2,2)) / (G(1,1)*G(2,2)-
     >G(1,2)*G(2,1)))
      F0(2,1) = F0(1,2)
      F0(3,4) = +(-(G(3,3)*G(4,3)+G(3,4)*G(4,4)) / (G(3,3)*G(4,4)-
     >G(3,4)*G(4,3)) )
      F0(4,3) = F0(3,4)
      F0(2,2) = (1.+F0(1,2)**2)/F0(1,1)
      F0(4,4) = (1.+F0(3,4)**2)/F0(3,3)

      RETURN       
      END
