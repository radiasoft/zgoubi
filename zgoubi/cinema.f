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
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE CINEMA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE "C.TRAJ.H"     ! COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI
      INCLUDE "C.CINE.H"     ! COMMON/CINE/ M1,M2,M3,M4,M12,M1212,M22,P0 ,G,C,C2,BG,EL3M,PC3,THETA,BETA,Q,PS(5),TS(5),NPS,NTS,II
      DOUBLE PRECISION M1,M2,M3,M4,M12,M1212,M22
 
      TF=THETA+TS(NTS)-BETA
      CS=SIN(P)*SIN(PS(NPS))+COS(TF-T)*COS(P)*COS(PS(NPS))
      CS2=CS*CS
      SS2=1.D0-CS2
      T=TF
      P=PS(NPS)
      BET=THETA-BETA
      CALL CHAREF(.FALSE.,ZERO,ZERO,BET)
      V=SS2+CS2/(G*G)
      DIS=SQRT(V-C2*SS2)/G
      U=(-C*SS2+CS*DIS)/V
      EL3=EL3M+BG*PC3*U
      PL3=SQRT(EL3*(EL3+2.D0*M3))
      DP=PL3/P0
      CALL MAJTRA(II)
      RETURN
      END
