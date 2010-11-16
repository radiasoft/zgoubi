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
      SUBROUTINE DRTENT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------
C     Droite de coupure entree
C-----------------------------
      COMMON/DROITE/ AM(9),BM(9),CM(9),IDRT
      COMMON/MARK/ KART,KALC,KERK,KUASEX
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG

      CT=COS(T)
      ST=SIN(T)
 
        DEN = AM(1)*CT+BM(1)*ST

        IF    (KART .EQ. 1) THEN
C         ... CARTESIEN
          XN = (-CM(1)*CT+X*BM(1)*ST-Y*BM(1)*CT)/DEN
          Y  = (-CM(1)*ST-X*AM(1)*ST+Y*AM(1)*CT)/DEN
          DL = (XN-X)/CT
          SAR = SAR + DL/COS(P)
          Z = Z + DL*TAN(P)
          X = XN
C          WRITE(NRES,*) ' SBR INTEG, DR. ENTREE , X,Y=',X,Y
        ELSEIF(KART .EQ. 2) THEN
C         ... POLAIRE
C ****************  NEVER BEEN USED! CHECK FORMULAE FIRST !!! *************8
          XN = (-CM(1)*ST+Y*AM(1)*CT)/DEN
          YN = (-CM(1)*CT-Y*BM(1)*CT)/DEN
          SAR = SAR + XN/CT*COS(P)
          Z = Z + XN/CT*TAN(P)
          X = ATAN(XN/YN)
          Y = YN/COS(X)
          T = T+X
        ENDIF
 
      RETURN
      END
