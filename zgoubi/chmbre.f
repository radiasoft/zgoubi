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
      SUBROUTINE CHMBRE(IT,Y,Z,SAR,
     >                             KEX,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------------------------------------
C        COMPTAGE DES PARTICULES QUI SORTENT D'UNE
C     CHAMBRE DE LIMITES TRANSVERSALES YLIM2 ET ZLIM2.
C     ------------------------------------------------
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YL2,ZL2,SORT(MXT),FMAG,BMAX
     > ,YC,ZC
 
C------- SKIP SI LA PARTICULE EST DEJA SORTIE :
      IF(KEX .LT. -1) RETURN
 
CC---- SEULEMNT SI LA FACE MAGNTQ EST FRANCHIE :
CCCCC IF(FMAG .GE. .45D0) THEN
 
        TEMP = SAR

        YP2 = (Y - YC)*(Y - YC)/YL2
        ZP2 = (Z - ZC)*(Z - ZC)/ZL2

        IF    (IFORM .EQ. 1) THEN
C--------- CHAMBRE RECTANGULAIRE
          IF( YP2 .GE. 1.D0 .OR. ZP2 .GE. 1.D0) CALL KSTOP(4,IT,KEX,*99)
        ELSEIF(IFORM .EQ. 2) THEN
C--------- CHAMBRE ELLIPTIQUE
          IF( ( YP2 + ZP2 ) .GE. 1.D0 ) CALL KSTOP(4,IT,KEX,*99)
        ENDIF

CCCCC ENDIF
 
      RETURN
 99   RETURN 1
      END
