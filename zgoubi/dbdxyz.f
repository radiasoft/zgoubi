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
      SUBROUTINE DBDXYZ(IDB,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DB(3,3),DDB(3,3,3)
      DIMENSION D3BX(3,3,3), D3BY(3,3,3), D3BZ(3,3,3)
      DIMENSION D4BX(3,3,3,3) ,D4BY(3,3,3,3) ,D4BZ(3,3,3,3)
C     -----------------------------------------
C     Compute  derivatives di+j+k(B)/dXidYjdZk.
C     B stands for magnetic or electric field.
C     -----------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
 
C----- DERIVEES 1-ERES DE BX, BY, BZ(X,Y,Z)
C     INDICE 1 = ELEMENTS DIFFERENTIELS : DX, DY, OU DZ.
C     INDICE 2 DU TABLEAU = COORDONNEE : BX, BY, OU BZ,
 
      DB(3,3)  = - ( DB(1,1) + DB(2,2) )
      DB(1,2) = DB(2,1)
      DB(1,3) = DB(3,1)
      DB(2,3) = DB(3,2)
 
      IF(IDB .EQ. 1) RETURN
 
C----- DERIVEES 2-EMES DE BX, BY, BZ(X,Y,Z).
C     INDICES 1,2 = ELEMENTS DIFFERENTIELS : DX, DY, OU DZ.
C     INDICE 3 DU TABLEAU = COORDONNEE : BX, BY, OU BZ,
 
      DDB(3,3,1) = -(DDB(1,1,1) + DDB(2,2,1))
      DDB(3,3,2) = -(DDB(2,1,1) + DDB(2,2,2))
      DDB(3,3,3) = -(DDB(3,1,1) + DDB(3,2,2))
 
      DDB(1,2,1) = DDB(2,1,1)
      DDB(1,3,1) = DDB(3,1,1)
      DDB(2,3,1) = DDB(3,2,1)
 
      DDB(1,1,2) = DDB(2,1,1)
      DDB(2,1,2) = DDB(2,2,1)
      DDB(3,1,2) = DDB(3,2,1)
      DDB(1,2,2) = DDB(2,2,1)
      DDB(1,3,2) = DDB(3,2,1)
      DDB(2,3,2) = DDB(3,2,2)
 
      DDB(1,1,3) = DDB(3,1,1)
      DDB(2,1,3) = DDB(3,2,1)
      DDB(3,1,3) = DDB(3,3,1)
      DDB(1,2,3) = DDB(3,2,1)
      DDB(2,2,3) = DDB(3,2,2)
      DDB(3,2,3) = DDB(3,3,2)
      DDB(1,3,3) = DDB(3,3,1)
      DDB(2,3,3) = DDB(3,3,2)
 
      IF(IDB .EQ. 2) RETURN
 
C----- 0 = dDivB/dXZ, dDivB/dYZ
      D3BX(3,3,3)   =  -(D3BX(3,1,1)+D3BX(3,2,2))
      D3BY(3,3,3)   =  -(D3BX(3,2,1)+D3BY(3,2,2))
 
C----- DERIVEES 3-EMES DE BX
      D3BX(1,3,1) = D3BX(3,1,1)
      D3BX(2,3,1) = D3BX(3,2,1)
      D3BX(3,1,2) = D3BX(3,2,1)
      D3BX(1,3,2) = D3BX(3,2,1)
      D3BX(2,3,2) = D3BX(3,2,2)
      D3BX(1,1,3) = D3BX(3,1,1)
      D3BX(2,1,3) = D3BX(3,2,1)
      D3BX(1,2,3) = D3BX(3,2,1)
      D3BX(2,2,3) = D3BX(3,2,2)
 
C     DERIVEES 3-EMES DE BY
      D3BY(3,1,1) = D3BX(3,2,1)
      D3BY(3,2,1) = D3BX(3,2,2)
      D3BY(1,3,1) = D3BX(3,2,1)
      D3BY(2,3,1) = D3BX(3,2,2)
      D3BY(3,1,2) = D3BX(3,2,2)
      D3BY(1,3,2) = D3BX(3,2,2)
      D3BY(2,3,2) = D3BY(3,2,2)
      D3BY(1,1,3) = D3BX(3,2,1)
      D3BY(2,1,3) = D3BX(3,2,2)
      D3BY(1,2,3) = D3BX(3,2,2)
      D3BY(2,2,3) = D3BY(3,2,2)
 
C     DERIVEES 3-EMES DE BZ
      D3BZ(1,1,1) = D3BX(3,1,1)
      D3BZ(2,1,1) = D3BX(3,2,1)
      D3BZ(1,2,1) = D3BX(3,2,1)
      D3BZ(2,2,1) = D3BX(3,2,2)
      D3BZ(3,3,1) = D3BX(3,3,3)
      D3BZ(1,1,2) = D3BX(3,2,1)
      D3BZ(2,1,2) = D3BX(3,2,2)
      D3BZ(1,2,2) = D3BX(3,2,2)
      D3BZ(2,2,2) = D3BY(3,2,2)
      D3BZ(3,3,2) = D3BY(3,3,3)
      D3BZ(3,1,3) = D3BX(3,3,3)
      D3BZ(3,2,3) = D3BY(3,3,3)
      D3BZ(1,3,3) = D3BX(3,3,3)
      D3BZ(2,3,3) = D3BY(3,3,3)
 
      IF(IDB .EQ. 3) RETURN
 
C----- 0 = dDivB/dX2, dDivB/dXY, dDivB/dY2, dDivB/dZ2
      D3BX(3,3,1)   =  -(D3BX(1,1,1)+D3BX(2,2,1))
      D3BX(3,3,2)   =  -(D3BX(2,1,1)+D3BX(2,2,2))
      D3BY(3,3,2)   =  -(D3BX(2,2,1)+D3BY(2,2,2))
      D3BZ(3,3,3)   =  -(D3BX(3,3,1)+D3BY(3,3,2))
 
C       DERIVEES 3-EMES DE BX
        D3BX(1,2,1) = D3BX(2,1,1)
        D3BX(1,1,2) = D3BX(2,1,1)
        D3BX(2,1,2) = D3BX(2,2,1)
        D3BX(1,2,2) = D3BX(2,2,1)
        D3BX(3,1,3) = D3BX(3,3,1)
        D3BX(3,2,3) = D3BX(3,3,2)
        D3BX(1,3,3) = D3BX(3,3,1)
        D3BX(2,3,3) = D3BX(3,3,2)
 
C       DERIVEES 3-EMES DE BY
        D3BY(1,1,1) = D3BX(2,1,1)
        D3BY(2,1,1) = D3BX(2,2,1)
        D3BY(1,2,1) = D3BX(2,2,1)
        D3BY(2,2,1) = D3BX(2,2,2)
        D3BY(3,3,1) = D3BX(3,3,2)
        D3BY(1,1,2) = D3BX(2,2,1)
        D3BY(2,1,2) = D3BX(2,2,2)
        D3BY(1,2,2) = D3BX(2,2,2)
        D3BY(3,1,3) = D3BX(3,3,2)
        D3BY(3,2,3) = D3BY(3,3,2)
        D3BY(1,3,3) = D3BX(3,3,2)
        D3BY(2,3,3) = D3BY(3,3,2)
 
C------- DERIVEES 3-EMES DE BZ
        D3BZ(3,1,1) = D3BX(3,3,1)
        D3BZ(3,2,1) = D3BX(3,3,2)
        D3BZ(1,3,1) = D3BX(3,3,1)
        D3BZ(2,3,1) = D3BX(3,3,2)
        D3BZ(3,1,2) = D3BX(3,3,2)
        D3BZ(3,2,2) = D3BY(3,3,2)
        D3BZ(1,3,2) = D3BX(3,3,2)
        D3BZ(2,3,2) = D3BY(3,3,2)
        D3BZ(1,1,3) = D3BX(3,3,1)
        D3BZ(2,1,3) = D3BX(3,3,2)
        D3BZ(1,2,3) = D3BX(3,3,2)
        D3BZ(2,2,3) = D3BY(3,3,2)
 
C      IF(IDB .EQ. 3) RETURN
C------- DERIVEES 4-EMES
 
        D4BZ(3,3,3,3)=D4BX(3,1,1,1) + 2.D0*D4BX(3,2,2,1) + D4BY(3,2,2,2)
 
C------- DERIVEES 4-EMES DE BX
        D4BX(1,3,1,1) = D4BX(3,1,1,1)
        D4BX(2,3,1,1) = D4BX(3,2,1,1)
 
        D4BX(3,1,2,1) = D4BX(3,2,1,1)
        D4BX(1,3,2,1) = D4BX(3,2,1,1)
        D4BX(2,3,2,1) = D4BX(3,2,2,1)
 
        D4BX(1,1,3,1) = D4BX(3,1,1,1)
        D4BX(2,1,3,1) = D4BX(3,2,1,1)
        D4BX(1,2,3,1) = D4BX(3,2,1,1)
        D4BX(2,2,3,1) = D4BX(3,2,2,1)
 
        D4BX(3,1,1,2) = D4BX(3,2,1,1)
        D4BX(3,2,1,2) = D4BX(3,2,2,1)
        D4BX(1,3,1,2) = D4BX(3,2,1,1)
        D4BX(2,3,1,2) = D4BX(3,2,2,1)
 
        D4BX(3,1,2,2) = D4BX(3,2,2,1)
        D4BX(1,3,2,2) = D4BX(3,2,2,1)
        D4BX(2,3,2,2) = D4BX(3,2,2,2)
 
        D4BX(1,1,3,2) = D4BX(3,2,1,1)
        D4BX(2,1,3,2) = D4BX(3,2,2,1)
        D4BX(1,2,3,2) = D4BX(3,2,2,1)
        D4BX(2,2,3,2) = D4BX(3,2,2,2)
 
        D4BX(1,1,1,3) = D4BX(3,1,1,1)
        D4BX(2,1,1,3) = D4BX(3,2,1,1)
        D4BX(1,2,1,3) = D4BX(3,2,1,1)
        D4BX(2,2,1,3) = D4BX(3,2,2,1)
        D4BX(3,3,1,3) = D4BX(3,3,3,1)
 
        D4BX(1,1,2,3) = D4BX(3,2,1,1)
        D4BX(2,1,2,3) = D4BX(3,2,2,1)
        D4BX(1,2,2,3) = D4BX(3,2,2,1)
        D4BX(2,2,2,3) = D4BX(3,2,2,2)
        D4BX(3,3,2,3) = D4BX(3,3,3,2)
 
        D4BX(3,1,3,3) = D4BX(3,3,3,1)
        D4BX(3,2,3,3) = D4BX(3,3,3,2)
        D4BX(1,3,3,3) = D4BX(3,3,3,1)
        D4BX(2,3,3,3) = D4BX(3,3,3,2)
 
C------- DERIVEES 4-EMES DE BY
        D4BY(3,1,1,1) = D4BX(3,2,1,1)
        D4BY(3,2,1,1) = D4BX(3,2,2,1)
        D4BY(1,3,1,1) = D4BX(3,2,1,1)
        D4BY(2,3,1,1) = D4BX(3,2,2,1)
 
        D4BY(3,1,2,1) = D4BX(3,2,2,1)
        D4BY(3,2,2,1) = D4BX(3,2,2,2)
        D4BY(1,3,2,1) = D4BX(3,2,2,1)
        D4BY(2,3,2,1) = D4BX(3,2,2,2)
 
        D4BY(1,1,3,1) = D4BX(3,2,1,1)
        D4BY(2,1,3,1) = D4BX(3,2,2,1)
        D4BY(1,2,3,1) = D4BX(3,2,2,1)
        D4BY(2,2,3,1) = D4BX(3,2,2,2)
        D4BY(3,3,3,1) = D4BX(3,3,3,2)
 
        D4BY(3,1,1,2) = D4BX(3,2,2,1)
        D4BY(3,2,1,2) = D4BX(3,2,2,2)
        D4BY(1,3,1,2) = D4BX(3,2,2,1)
        D4BY(2,3,1,2) = D4BX(3,2,2,2)
 
        D4BY(3,1,2,2) = D4BX(3,2,2,2)
        D4BY(1,3,2,2) = D4BX(3,2,2,2)
        D4BY(2,3,2,2) = D4BY(3,2,2,2)
 
        D4BY(1,1,3,2) = D4BX(3,2,2,1)
        D4BY(2,1,3,2) = D4BX(3,2,2,2)
        D4BY(1,2,3,2) = D4BX(3,2,2,2)
        D4BY(2,2,3,2) = D4BY(3,2,2,2)
 
        D4BY(1,1,1,3) = D4BX(3,2,1,1)
        D4BY(2,1,1,3) = D4BX(3,2,2,1)
        D4BY(1,2,1,3) = D4BX(3,2,2,1)
        D4BY(2,2,1,3) = D4BX(3,2,2,2)
        D4BY(3,3,1,3) = D4BX(3,3,3,2)
 
        D4BY(1,1,2,3) = D4BX(3,2,2,1)
        D4BY(2,1,2,3) = D4BX(3,2,2,2)
        D4BY(1,2,2,3) = D4BX(3,2,2,2)
        D4BY(2,2,2,3) = D4BY(3,2,2,2)
        D4BY(3,3,2,3) = D4BY(3,3,3,2)
 
        D4BY(3,1,3,3) = D4BX(3,3,3,2)
        D4BY(3,2,3,3) = D4BY(3,3,3,2)
        D4BY(1,3,3,3) = D4BX(3,3,3,2)
        D4BY(2,3,3,3) = D4BY(3,3,3,2)
 
C------- DERIVEES 4-EMES DE BZ
        D4BZ(1,1,1,1) = D4BX(3,1,1,1)
        D4BZ(2,1,1,1) = D4BX(3,2,1,1)
        D4BZ(1,2,1,1) = D4BX(3,2,1,1)
        D4BZ(2,2,1,1) = D4BX(3,2,2,1)
        D4BZ(3,3,1,1) = D4BX(3,3,3,1)
 
        D4BZ(1,1,2,1) = D4BX(3,2,1,1)
        D4BZ(2,1,2,1) = D4BX(3,2,2,1)
        D4BZ(1,2,2,1) = D4BX(3,2,2,1)
        D4BZ(2,2,2,1) = D4BX(3,2,2,2)
        D4BZ(3,3,2,1) = D4BX(3,3,3,2)
 
        D4BZ(3,1,3,1) = D4BX(3,3,3,1)
        D4BZ(3,2,3,1) = D4BX(3,3,3,2)
        D4BZ(1,3,3,1) = D4BX(3,3,3,1)
        D4BZ(2,3,3,1) = D4BX(3,3,3,2)
 
        D4BZ(1,1,1,2) = D4BX(3,2,1,1)
        D4BZ(2,1,1,2) = D4BX(3,2,2,1)
        D4BZ(1,2,1,2) = D4BX(3,2,2,1)
        D4BZ(2,2,1,2) = D4BX(3,2,2,2)
        D4BZ(3,3,1,2) = D4BX(3,3,3,2)
 
        D4BZ(1,1,2,2) = D4BX(3,2,2,1)
        D4BZ(2,1,2,2) = D4BX(3,2,2,2)
        D4BZ(1,2,2,2) = D4BX(3,2,2,2)
        D4BZ(2,2,2,2) = D4BY(3,2,2,2)
        D4BZ(3,3,2,2) = D4BY(3,3,3,2)
 
        D4BZ(3,1,3,2) = D4BX(3,3,3,2)
        D4BZ(3,2,3,2) = D4BY(3,3,3,2)
        D4BZ(1,3,3,2) = D4BX(3,3,3,2)
        D4BZ(2,3,3,2) = D4BY(3,3,3,2)
 
        D4BZ(3,1,1,3) = D4BX(3,3,3,1)
        D4BZ(3,2,1,3) = D4BX(3,3,3,2)
        D4BZ(1,3,1,3) = D4BX(3,3,3,1)
        D4BZ(2,3,1,3) = D4BX(3,3,3,2)
 
        D4BZ(3,1,2,3) = D4BX(3,3,3,2)
        D4BZ(3,2,2,3) = D4BY(3,3,3,2)
        D4BZ(1,3,2,3) = D4BX(3,3,3,2)
        D4BZ(2,3,2,3) = D4BY(3,3,3,2)
 
        D4BZ(1,1,3,3) = D4BX(3,3,3,1)
        D4BZ(2,1,3,3) = D4BX(3,3,3,2)
        D4BZ(1,2,3,3) = D4BX(3,3,3,2)
        D4BZ(2,2,3,3) = D4BY(3,3,3,2)
 
      RETURN
      END
