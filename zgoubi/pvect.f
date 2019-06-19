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
C  -------
      SUBROUTINE PVECT(I,J,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ----------------------------------------------
C     CALCUL DE LA DERIVEE DU PRODUIT VECTORIEL UxB
C     LE 1-ER  INDICE EST L'ORDRE DE DERIVATION /DS,
C     LE 2-EME INDICE EST LE NUMERO DE COORDONNEE,
C     EXEMPLE: B(2,3)=dBz/ds, V(3,2)=d2Vy/ds2
C     ----------------------------------------------
      INCLUDE "C.VITES.H"     ! COMMON/VITES/ U(6,3),DQBR(6),DDT(6)
      INCLUDE "C.CHAVE.H"     ! COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
 
      V(I,1)=U(J,2)*B(K,3)-U(J,3)*B(K,2)
      V(I,2)=U(J,3)*B(K,1)-U(J,1)*B(K,3)
      V(I,3)=U(J,1)*B(K,2)-U(J,2)*B(K,1)
 
      RETURN
      END
