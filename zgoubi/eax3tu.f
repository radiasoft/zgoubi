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
      SUBROUTINE EAX3TU(OR,X,VDR,D,X22,EX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EX(*)
 
C------- LENTILLE A 3 TUBES DISTANTS DE D.
C        V(X) = (V2-V1)/2.Omega.D ( LOG(CH++/CH+) + LOG(CH--/CH-) )
 
        EPP=EXP(OR*(X+X22+D))
        EP=EXP(OR*(X+X22))
        EMM=EXP(OR*(X-X22-D))
        EM=EXP(OR*(X-X22))

        CPP=.5D0*(EPP+1.D0/EPP)
        CP=.5D0*(EP+1.D0/EP)
        CMM=.5D0*(EMM+1.D0/EMM)
        CM=.5D0*(EM+1.D0/EM)
        CPP2=CPP*CPP
        CP2=CP*CP
        CPP4=CPP2*CPP2
        CP4=CP2*CP2
        CMM2=CMM*CMM
        CM2=CM*CM
        CMM4=CMM2*CMM2
        CM4=CM2*CM2

        SPP=.5D0*(EPP-1.D0/EPP)
        SP=.5D0*(EP-1.D0/EP)
        SMM=.5D0*(EMM-1.D0/EMM)
        SM=.5D0*(EM-1.D0/EM)
        SPP2=SPP*SPP
        SP2=SP*SP
        SPP4=SPP2*SPP2
        SP4=SP2*SP2
        SMM2=SMM*SMM
        SM2=SM*SM
        SMM4=SMM2*SMM2
        SM4=SM2*SM2
 
        OR2=OR*OR

C------- EX(n)=-dnV/dXn
        EX(1)=-.5D0*VDR*(SPP/CPP - SP/CP + SMM/CMM - SM/CM) 
        EX(2)=-OR*.5D0*VDR*(1.D0/CPP2 - 1.D0/CP2 + 1.D0/CMM2 - 1.D0/CM2)
        EX(3)= OR2*VDR*(SPP/CPP2/CPP - SP/CP2/CP 
     >     + SMM/CMM2/CMM - SM/CM2/CM) 
        EX(4)= OR2*OR*VDR*( (1.D0-2.D0*SPP2)/CPP4 - (1.D0-2.D0*SP2)/CP4
     >   + (1.D0-2.D0*SMM2)/CMM4 - (1.D0-2.D0*SM2)/CM4  )
        EX(5)=-4.D0*OR2*OR2*VDR*( SPP*(2.D0-SPP2)/CPP4/CPP - 
     >     SP*(2.D0-SP2)/CP4/CP + SMM*(2.D0-SMM2)/CMM4/CMM -  
     >     SM*(2.D0-SM2)/CM4/CM )
C---------- EX(6) still waiting for provision...
        EX(6)=0D0

      RETURN
      END
