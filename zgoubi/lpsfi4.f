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
      SUBROUTINE LPSFI4(
     >                  SQX,SQZ,S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(4,4)

      INCLUDE "C.CONST_3.H"      ! COMMON/CONST/ CL9,CEL,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
C----- CONVERSION DES COORD. (CM,MRD) -> (M,RD)
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)

        UNIT1=  UNIT(1)
        UNIT2=  UNIT(2)
        UNIT3=  UNIT(3)
        UNIT4=  UNIT(4)

        XM=0.D0
        XPM=0.D0
        ZM=0.D0
        ZPM=0.D0
        NPTS = 0
        DO 21 I=1,IMAX
          IF(IEX(I) .LT. -1) GOTO 21
          NPTS = NPTS + 1
          X   = F(2,I)*UNIT1
          XP  = F(3,I)*UNIT2
          Z   = F(4,I)*UNIT3
          ZP  = F(5,I)*UNIT4
          XM  = XM  + X
          XPM = XPM + XP
          ZM  = ZM  + Z
          ZPM = ZPM + ZP
 21     CONTINUE
        XM =  XM/NPTS
        XPM = XPM/NPTS
        ZM =  ZM/NPTS
        ZPM = ZPM/NPTS
        UM= XM
        UPM=XPM
        VM= ZM
        VPM=ZPM

        X2  = 0.D0
        XP2 = 0.D0
        XXP = 0.D0
        Z2  = 0.D0
        ZP2 = 0.D0
        ZZP = 0.D0
        XZ  = 0.D0
        XZP = 0.D0
        XPZ = 0.D0
        XPZP= 0.D0
        DO 26 I=1,IMAX
          IF(IEX(I) .LT. -1) GOTO 26
          X   = F(2,I)*UNIT1
          XP  = F(3,I)*UNIT2
          Z   = F(4,I)*UNIT3
          ZP  = F(5,I)*UNIT4
          X2  = X2  + (X- UM)**2
          XP2 = XP2 + (XP-UPM)**2
          XXP = XXP + (X- UM)*(XP-UPM)
          Z2  = Z2  + (Z- VM)**2
          ZP2 = ZP2 + (ZP-VPM)**2
          ZZP = ZZP + (Z- VM)*(ZP-VPM)
          XZ  = XZ  + (X- UM)*(Z-VM)
          XZP = XZP + (X- UM)*(ZP-VPM)
          XPZ = XPZ + (XP- UPM)*(Z-VM)
          XPZP = XPZP + (XP- UPM)*(ZP-VPM)
 26     CONTINUE
        S(1,1) = X2/NPTS
        S(1,2) = XXP/NPTS
        S(1,3) = XZ /NPTS
        S(1,4) = XZP/NPTS
        S(2,1) = S(1,2)
        S(2,2) = XP2/NPTS
        S(2,3) = XPZ/NPTS
        S(2,4) = XPZP/NPTS
        S(3,1) = S(1,3)
        S(3,2) = S(2,3)
        S(3,3) = Z2/NPTS
        S(3,4) = ZZP/NPTS
        S(4,1) = S(1,4)
        S(4,2) = S(2,4)
        S(4,3) = S(3,4)
        S(4,4) = ZP2/NPTS

C G. Leleux : surface de l'ellipse S=4.pi.sqrt(DELTA)
C Soit d11=X2/sqrt(DELTA), d12=XXP/sqrt(DELTA), d22=XP2/sqrt(DELTA), alors 
C d22.x^2-2.d12.x.x'+d11.x'^2=S/pi=4sqrt(DELTA), ce qui permet d'ecrire 
C gamma=d22=XP2/sqrt(DELTA),-alpha=d12=XXP/sqrt(DELTA),beta=d11=X2/sqrt(DELTA).
C En outre, par definition des dij, 
C     2.sigma_x=sqrt(d11.S/pi),  2.sigma_x'=sqrt(d22.S/pi). 
C En outre, frontiere : 
C          <x^2>_frontiere=2.(sigma_x)^2,    <x'^2>_frontiere=2.(sigma_x')^2

C------- Courant invariant at 1 sigma is U=4.sqrt(DELTA)=Eps/pi (consistant with zgoubi !!) :
C Eps=ellipse surface
        SQX = SQRT(S(1,1)*S(2,2)-S(1,2)*S(1,2)) 
        SQZ = SQRT(S(3,3)*S(4,4)-S(3,4)*S(3,4)) 

      RETURN
      END
