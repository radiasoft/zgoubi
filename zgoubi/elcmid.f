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
C  Brookhaven National Laboratory               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE ELCMID(R,BRI,
     >                           ER0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C--------------------------------------------------------------
C  Compute ER at Z=0 & derivatives wrt to R, as a function of R
C--------------------------------------------------------------
      DIMENSION ER0(*)

      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M

      SAVE D, D2, D3, D4, R1, R2, V21D, V32D

      V21 = V21D * BRI
      V32 = V32D * BRI
      PIV21 = PI * V21
      PIV32 = PI * V32 
      PI2V21 = PI*PIV21
      PI2V32 = PI*PIV32
      PI3V21 = PI*PI2V21
      PI3V32 = PI*PI2V32
      PI4V21 = PI*PI3V21
      PI4V32 = PI*PI3V32

C------- POUR VERIF CONSERVATION DE L'ENERGIE AU POINT COURANT
C      VR =  0.5D0 * (V21D + V32D)*D + 
C     > ((V21D*D)*ATAN(SINH((PI*(R - R1))/D)))/PI + 
C     >  ((V32D*D)*ATAN(SINH((PI*(R - R2))/D)))/PI
C         IF(R .LT. 24) WRITE(6,*) R,VR
C      CALL ENRGW(R,VR)
C--------------------------------

C------- ER0(N)=-DNV/DRN
C------ ER0 : 
          SEC1 = 1.D0/COSH((PI*(R - R1))/D)
          SEC2 = 1.D0/COSH((PI*(R - R2))/D)
      ER0(1) =  -(V21*SEC1 + V32*SEC2)

C------ DER0/DR :
          TAN1 = TANH((PI*(R - R1))/D)
          TAN2 = TANH((PI*(R - R2))/D)
      ER0(2) =   (PIV21*SEC1*TAN1 + PIV32*SEC2*TAN2)/D
C------ D2ER0/DR2 :
      SEC12 = SEC1 * SEC1
      SEC22 = SEC2 * SEC2
      ER0(3) = (PI2V21*SEC1*(-1.D0 + 2.D0*SEC12) +  
     >          PI2V32*SEC2*(-1.D0 + 2.D0*SEC22)) /D2
C------ D3ER0/DR3 :
      ER0(4) = (PI3V21*SEC1*(1.D0 - 6.D0*SEC12)*TAN1 + 
     >          PI3V32*SEC2*(1.D0 - 6.D0*SEC22)*TAN2)/ D3
C------ D4ER0/DR4 :
      ER0(5) = 
     >((-(PI4V21*(115 - 76*COSH((2.D0*PI*(R - R1))/D) + 
     >                      COSH((4.D0*PI*(R - R1))/D))*SEC1**5) - 
     >    PI4V32*(115 - 76*COSH((2.D0*PI*(R - R2))/D) + 
     >                      COSH((4.D0*PI*(R - R2))/D))*SEC2**5
     >      ))/(8.D0*D4)

      RETURN

      ENTRY ELCMIW(DIN,R1IN,R2IN,V21DIN,V32DIN)
      D = DIN
      D2 = D*D
      D3 = D*D2
      D4 = D2*D2
      R1 = R1IN 
      R2 = R2IN 
      V21D = V21DIN 
      V32D = V32DIN 

      RETURN
      END
