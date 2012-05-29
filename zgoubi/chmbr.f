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
      SUBROUTINE CHMBR(I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------------------------------------
C        COMPTAGE DES PARTICULES QUI SORTENT D'UNE
C     CHAMBRE DE LIMITES TRANSVERSALES YLIM2 ET ZLIM2.
C     ------------------------------------------------
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
 
      IT = I1-1
 2    CONTINUE
        IT=IT+1
        IF( IEX(IT) .LT. -1) GOTO 3
        CALL CHMBRE(IT,F(2,IT),F(4,IT),F(6,IT),
     >                                         IEX(IT),*3)
 3      CONTINUE
      IF(IT .LT. I2) GOTO 2
 
      RETURN
      END
