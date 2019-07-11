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
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE EAXIAL(OM,OX,V21,D,
     >                              EX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EX(*)

      IF(D .EQ. 0.D0) THEN
C------- LENTILLE A 2 TUBES INFINIMENT PROCHES, cf SEPTIER
C        V(X) = .5(V1+V2) ( 1 + (V2-V1)/(V2+V1)TH(OmegaX) )

        E=EXP(OX)
        CH=.5D0*(E+1.D0/E)
        SH=.5D0*(E-1.D0/E)
        TH=SH/CH

C------- EX(n)=-dnV/dXn
        EX(1)=-.5D0*V21*OM/CH/CH
        EX(2)=-OM* 2.D0*TH*EX(1)
        EX(3)=-OM*(3.D0*TH*EX(2)+  2.D0*EX(1)*OM )
        EX(4)=-OM*(4.D0*TH*EX(3)+( 5.D0*EX(2)+   2.D0*TH*EX(1)*OM )*OM)
        EX(5)=-OM*(5.D0*TH*EX(4)+( 9.D0*EX(3)+ ( 7.D0*TH*EX(2)
     >                   +  2.D0*EX(1)*OM )*OM)*OM )
        EX(6)=-OM*(6.D0*TH*EX(5)+(14.D0*EX(4)+ (16.D0*TH*EX(3)
     >                   + (9.D0*EX(2)+ 2.D0*TH*EX(1)*OM )*OM)*OM)*OM )

      ELSE
C------- LENTILLE A 2 TUBES DISTANTS DE D, cf SEPTIER
C        V(X) = .5(V1+V2) ( 1 + (V2-V1)/(V2+V1)/(2 OM D)*LOG(CH+/CH-) )

        EP=EXP(OX+OM*D)
        EM=EXP(OX-OM*D)
        CP=.5D0*(EP+1.D0/EP)
        CM=.5D0*(EM+1.D0/EM)
        CP2=CP*CP
        CM2=CM*CM
        SP=.5D0*(EP-1.D0/EP)
        SM=.5D0*(EM-1.D0/EM)
        SP2=SP*SP
        SM2=SM*SM
        TP=SP/CP
        TM=SM/CM
        VD=.25D0*V21/D

        EX(1)=  -VD*(TP-TM)
        VD=VD*OM
        EX(2)= -VD*(1.D0/CP2-1.D0/CM2)
        EX(3)= -2.D0*OM*( (TP+TM)*EX(2) + VD*(TM/CP2-TP/CM2) )
        EX(4)= -2.D0*OM*( OM*(1.D0/CP2+1.D0/CM2)*EX(2) + (TP+TM)*EX(3)
     >         +2.D0*VD*OM*TM*TP*( 1.D0/CM2-1.D0/CP2)  )
        EX(5)= 4.D0*OM*OM*(  OM*(TM/CM2+TP/CP2)*EX(2)
     >  - (1.D0/CM2 + 1.D0/CP2)*EX(3) - (TP+TM)*EX(4)/2.D0/OM
     >  + VD*OM*( (1.D0/CM2-1.D0/CP2)*(TM/CP2+TP/CM2)
     >            - 2.D0*TP*TM*(TM/CM2-TP/CP2)) )
C-------- EX(6) waiting to be documented...
        EX(6)=0D0
      ENDIF

      RETURN
      END
