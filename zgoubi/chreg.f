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
      FUNCTION CHREG(KART,X,Y,XE,AREG,BREG,CREG,DXI,TPAS,
     >                                                   KREG,AL,BL,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CHREG
C------- Exist entrance and/or exit fringe field regions and,
C        XPAS is coded as x.yyyE10 (constant  step) or as 
C                                         x.yyyV1aa (variable step)
C------------------------------------------------------------------
      DIMENSION AREG(*),BREG(*),CREG(*),TPAS(*)

        IF(KREG .LE. 2) THEN
          KRG0 = KREG
          IF( KART .EQ. 1) THEN
C--------- COORDONNEES CARTESIENNES
            IF(KRG0 .EQ. 1) THEN
C----------- Now in entrance fringe field
C FM 08/99              D = 2.D0*XE - X
              D = -(AREG(1)*X + BREG(1)*Y + CREG(1))
C----------- Test End of entrance fringe-field:
              IF( D.LT.DXI ) KREG = 2
            ELSEIF(KRG0 .EQ. 2) THEN
C------------- Now in central region
C FM 08/99              D = 2.D0*XS-XLIM - X
              D = -(AREG(2)*X + BREG(2)*Y + CREG(2))
C		........................1            
C(1) Warning: Variable Y is used before its value has been defined
              IF( D.LE.DXI ) THEN
C--------------- End of central region
                IF(TPAS(3) .EQ. 0.D0) THEN
C----------------- Sharp edge at exit
                  KREG = 2
                ELSE
C----------------- Exit fringe field
                  KREG=3
                ENDIF
              ENDIF
            ENDIF
          ENDIF
          CHREG = KREG.NE.KRG0
          IF(CHREG) THEN
            AL = AREG(KRG0)
            BL = BREG(KRG0)
          ENDIF
        ELSE
C--------- There are only 3 possible regions : entrance|body|exit !
          CHREG = .FALSE.
        ENDIF
      RETURN
      END
