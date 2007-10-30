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
      SUBROUTINE SYMMED(Z,IDZ,BZ0,
     >                           B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BZ0(5,5)
      DIMENSION B(5,3),DB(3,3),DDB(3,3,3)
      DIMENSION D3BX(3,3,3), D3BY(3,3,3), D3BZ(3,3,3)
      DIMENSION D4BX(3,3,3,3) ,D4BY(3,3,3,3) ,D4BZ(3,3,3,3)
C     ---------------------------------------------------
C     CARTES  2-D  OU  ELEMENTS  DEFINIS  PAR  BZ(X,Y,0):
C     SYMMETRISE PAR RAPPORT AU PLAN MEDIAN
C     ---------------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG

C Reminicences from ancient times
      B0=BZ0(1,1) 
      B1=BZ0(2,1)
      B2=BZ0(1,2)
      B11=BZ0(3,1)
      B12=BZ0(2,2)
      B22=BZ0(1,3)
      B111=BZ0(4,1)
      B112=BZ0(3,2)
      B122=BZ0(2,3)
      B222=BZ0(1,4)
      B1111=BZ0(5,1)
      B1112=BZ0(4,2)
      B1122=BZ0(3,3)
      B1222=BZ0(2,4)
      B2222=BZ0(1,5)

      B33=-(B11+B22)
 
CALCUL DES DERIVEES DE BX, BY, BZ(X,Y,Z) PAR DEVELOPPEMENTS DE TAYLOR
C EN Z. LES DERIVEES PAR RAPPORT A Z SONT EXPRIMEES EN TERMES DE
C COMBINAISONS DES DERIVEES PAR RAPPORT A X ET Y DU PLAN MEDIAN.
C      Z=Z1
      ZB33 = Z*B33
      DDB(1,1,1) = Z * B111
      DDB(2,1,1) = Z * B112
      DDB(2,2,1) = Z * B122
      DDB(2,2,2) = Z * B222
 
      ZZB12 = -Z * (DDB(1,1,1) + DDB(2,2,1))
      ZZB21 = -Z * (DDB(2,1,1) + DDB(2,2,2))
      DB(1,1) = Z*B11
      DB(2,1) = Z*B12
      DB(2,2) = Z*B22
      DB(3,1) = B1 + 0.5D0 * ZZB12
      DB(3,2) = B2 + 0.5D0 * ZZB21
 
      DDB(3,1,1)=B11
      DDB(3,2,1)=B12
      DDB(3,2,2)=B22

C----- CALCUL DE BX, BY, BZ(X,Y,Z) PAR EXTRAPOLATION EN Z
      B(1,1) = Z * B1 + Z * ZZB12 /6.D0
      B(1,2) = Z * B2 + Z * ZZB21 /6.D0
      B(1,3) = B0 + Z *  ZB33 * .5D0
 
      IF(IDZ.GE.3) THEN
C--------- EXTRAPOLATION A L'ORDRE 3 EN Z
        D3BX(3,1,1) = B111
        D3BX(3,2,1) = B112
        D3BX(3,2,2) = B122
        D3BY(3,2,2) = B222

        IF(IDZ.GE.4) THEN
 
C--------- EXTRAPOLATION A L'ORDRE 4 EN Z
          D4BX(3,1,1,1) = B1111
          D4BX(3,2,1,1) = B1112
          D4BX(3,2,2,1) = B1122
          D4BX(3,2,2,2) = B1222
          D4BY(3,2,2,2) = B2222
 
          D4BX(3,3,3,1)=-(B1111+B1122)
          D4BX(3,3,3,2)=-(B1112+B1222)
          D4BY(3,3,3,2)=-(B1122+B2222)
 
          D3BX(2,1,1)= Z*B1112
          D3BX(2,2,1)= Z*B1122
          D3BX(2,2,2)= Z*B1222
          D3BZ(3,3,3)=Z*(B1111+2.D0*B1122+B2222)
          D3BX(3,3,1)=Z*D4BX(3,3,3,1)
          D3BX(3,3,2)=Z*D4BX(3,3,3,2)
          D3BY(3,3,2)=Z*D4BY(3,3,3,2)
 
          Z2B422=Z*D3BX(3,3,1)*.5D0
          Z2B313=Z*D3BX(3,3,2)*.5D0
          Z2B224=Z*D3BY(3,3,2)*.5D0
 
          Z3B44=Z*Z*D3BZ(3,3,3)
          Z3B422=Z*Z2B422/6.D0
          Z3B313=Z*Z2B313/6.D0
          Z3B224=Z*Z2B224/6.D0
 
          B(1,3) = B(1,3) + Z * Z3B44 /24.D0
          DB(1,1) = DB(1,1) + Z3B422
          DB(2,1) = DB(2,1) + Z3B313
          DB(2,2) = DB(2,2) + Z3B224
 
          DDB(3,1,1)=B11+Z2B422
          DDB(3,2,1)=B12+Z2B313
          DDB(3,2,2)=B22+Z2B224
 
          D3BX(1,1,1)=Z*B1111
          D3BY(2,2,2)=Z*B2222
        ENDIF
      ENDIF

      RETURN
      END
