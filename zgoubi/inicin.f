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
      SUBROUTINE INICIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/RIGID/ BORO,DPREF,DP,BR
      COMMON/CINE/M1,M2,M3,M4,M12,M1212,M22,P0 ,G,C,C2,BG,EL3M,PC3,THETA
     U,BETA,Q,PS(5),TS(5),NPS,NTS,II
      DOUBLE PRECISION M1,M2,M3,M4,M12,M1212,M22
      PL2=DP*P0
      PL22=PL2*PL2
      EL2=PL22/(M2+SQRT(M22+PL22))
      WLT=M12+EL2
      WCT=SQRT(M1212+2.D0*EL2*M1)
      G=WLT/WCT
      B=PL2/WLT
      ECI=2.D0*EL2*M1/(WCT+M12)
      Q=Q*.001D0
      M4=M12-M3-Q
      ECF=ECI+Q
      EC3=ECF*(ECF+2.D0*M4)/(2.D0*WCT)
      PC3=SQRT(EC3*(EC3+2.D0*M3))
      BC3=PC3/(EC3+M3)
      C=B/BC3
      C2=C*C
      BG=B*G
      F=BG*BG/(1.D0+G)
      EL3M=EC3+F*(EC3+M3)
      RETURN
      END
