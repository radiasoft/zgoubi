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
      SUBROUTINE SREF(KSTP,FE,ITEM,
     >                             EVB,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ----------------------
C     Calculate SR  E-field 
C     ----------------------
      DIMENSION EVB(*)
C----- E is the radiated E-field in observer time

      COMMON/CONST/ CL9,CL,PI,RAD,DEG,QE,AMPROT, CM2M
      INCLUDE 'MXSTEP.H'
      INCLUDE 'MAXTRA.H'
      INCLUDE 'CSR.H'
C      COMMON/CSR/ KTRA,KCSR,YZXB(MXSTEP,41,36),DWC(MXT)
      INCLUDE "MAXCOO.H"
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AM ,Q 
      COMMON/UNITS/ UNIT(MXJ)

C-----          Emitor Emitor-Obs.                       Nrmlzd
C               posit.  direct     beta  beta-dot  field obs. pos. 
      DIMENSION RP(3),  RO(3),     BE(3), BED(3),   B(3), ON(3),  OX(3)
      DIMENSION E(3), ONO(3), BEO(3)

      SAVE RP, RO, OX, RR

C      DIMENSION ONB(3)

        B(1) = YZXB(KSTP,ITEM,30)                ! Bx, Tesla
        B(2) = YZXB(KSTP,ITEM,31)                ! By
        B(3) = YZXB(KSTP,ITEM,32)                ! Bz

        CT = COS( YZXB(KSTP,ITEM,3) )
        CP = COS( YZXB(KSTP,ITEM,5) )

C------- Normalized ( n(t') )
        ON(1) = RO(1)/RR
        ON(2) = RO(2)/RR
        ON(3) = RO(3)/RR

C------- Beta
C   yzxb(1)=-1+dp/p. Rigidity (T.m),  Momentum (MeV/c) 
        BRO = YZXB(KSTP,ITEM,36) * (1.D0 + YZXB(KSTP,ITEM,1)) 
        P0 = BRO * CL * 1.D-6 * Q    
        BTA = P0 / SQRT( P0*P0 + AM*AM )  
C   q/m.gamma, MKSA units
        QMG = BTA * CL * CL / (P0 * 1.D6)* Q 

        BE(1) = BTA * CT * CP
        BE(2) = BTA * SIN( YZXB(KSTP,ITEM,3) ) * CP
        BE(3) = BTA * SIN( YZXB(KSTP,ITEM,5) )
                 
C------- Acceleration = Beta-dot = q/m beta x Field
        BED(1) = QMG * ( BE(2)*B(3) - BE(3)*B(2) )
        BED(2) = QMG * ( BE(3)*B(1) - BE(1)*B(3) )
        BED(3) = QMG * ( BE(1)*B(2) - BE(2)*B(1) )

C------- 1 - n.Beta
        UNB = (1.D0 - ON(1)*BE(1)) - (ON(2)*BE(2)+ON(3)*BE(3) )

        G1 = BTA * AM / P0
          UNB2 = UNB * UNB
          UNB3 = UNB2 * UNB
C--------- (n - Beta) x Beta-dot
          V1 = (ON(2) - BE(2)) * BED(3) - (ON(3) - BE(3)) * BED(2)
          V2 = (ON(3) - BE(3)) * BED(1) - (ON(1) - BE(1)) * BED(3)
          V3 = (ON(1) - BE(1)) * BED(2) - (ON(2) - BE(2)) * BED(1)

C------- Radiated electric field x,y,z components ( V/m )
        E(1) = FE * ( ON(2)*V3 - ON(3)*V2 ) / UNB3 / RR
        E(2) = FE * ( ON(3)*V1 - ON(1)*V3 ) / UNB3 / RR
        E(3) = FE * ( ON(1)*V2 - ON(2)*V1 ) / UNB3 / RR
C------- EVB = E + vXB
        EVB(1) = ( ON(2)*E(3) - ON(3)*E(2) )
        EVB(2) = ( ON(3)*E(1) - ON(1)*E(3) )
        EVB(3) = ( ON(1)*E(2) - ON(2)*E(1) )
        EVB(1) = E(1) + ( BE(2)*EVB(3) - BE(3)*EVB(2) )
        EVB(2) = E(2) + ( BE(3)*EVB(1) - BE(1)*EVB(3) )
        EVB(3) = E(3) + ( BE(1)*EVB(2) - BE(2)*EVB(1) )

      RETURN      

      ENTRY SREF1(KSTP,ITE,
     >                     RRR,ONO,BEO,BTAO)
C----- Position of front particle (the observer, victim of the interaction)
        DS = FO(6,IT)                 ! cm  
        DX = DS * COS(T) * COS(P)  
        DY = DS * SIN(T) * COS(P)  
        DZ = DS * SIN(P)           
        OX(1) = (X + DX) * UNIT(5)    ! m
        OX(2) = (Y + DY) * UNIT(1) 
        OX(3) = (Z + DZ) * UNIT(3)  
C------- Position of emitor (tag=ITE) in magnet frame  ( r(t') )
        DS = FO(6,ITE) * UNIT(5)      ! m  
        CT = COS( YZXB(KSTP,ITE,3) )
        CP = COS( YZXB(KSTP,ITE,5) )
        ST = SIN( YZXB(KSTP,ITE,3) )
        SP = SIN( YZXB(KSTP,ITE,5) )
        DX = DS * CT * CP             ! m  
        DY = DS * ST * CP
        DZ = DS * SP
        RP(1) = YZXB(KSTP,ITE,7) + DX                 ! X (m)
        RP(2) = YZXB(KSTP,ITE,2) + DY                 ! Y (m)
        RP(3) = YZXB(KSTP,ITE,4) + DZ                 ! Z (m)
C------- Emitor to observer vector  ( R(t') )
        RO(1) = OX(1) - RP(1)
        RO(2) = OX(2) - RP(2)
        RO(3) = OX(3) - RP(3)
        RRR = SQRT(RO(1)*RO(1) + RO(2)*RO(2) + RO(3)*RO(3))
        RR = RRR
C------- Normalized ( n(t') )
        ONO(1) = RO(1)/RRR
        ONO(2) = RO(2)/RRR
        ONO(3) = RO(3)/RRR
C------- Beta emittor
C   yzxb(1)=dp/p. Rigidity (T.m),  Momentum (MeV/c) 
        BRO = YZXB(KSTP,ITE,36) * (1.D0 + YZXB(KSTP,ITE,1)) 
        P0 = BRO * CL * 1.D-6 * Q    
        BTA = P0 / SQRT( P0*P0 + AM*AM )  
        BTAO = BTA
        BEO(1) = BTA * CT * CP
        BEO(2) = BTA * ST * CP
        BEO(3) = BTA * SP
                 
      RETURN
      END
