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
C  Brookhaven National Laboratory                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE SRPLT(NLOG,KEP,DCX,CX,CY,IX,IY,UNB,RP,TPART,
     >                                     NOC,NPTS,NRMA,NRMC,
     >                                         TM,YM,TC,SYDX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C------- Compute scales (KEP=1), or Plot (KEP=2)
      DIMENSION RP(*)
      DIMENSION CX(*), CY(*)

      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      LOGICAL INRANG
      SAVE X

      LOGICAL NONEWT
      SAVE NONEWT
      DATA NONEWT / .FALSE. /

      IF(NONEWT) THEN
        X = X + DCX
      ELSE
        X = DCX
        NONEWT = .TRUE. 
      ENDIF

      IF(NOC .GT. 0) THEN

        IF( IY .EQ. 11 ) THEN   
          CY(IY) = UNB                      ! 1 - n.beta
        ELSEIF( IY .EQ. 10 ) THEN   
          CY(IY) = TPART                      ! t'
        ELSEIF( IY .GE. 7 ) THEN
          CY(IY) = RP(IY-6)                 ! Particle X, Y, Z position (m)
        ENDIF

        IF( KEP .EQ. 1 ) THEN
C----------- Calculate scales

          IF(NOC .GE. 2) THEN

            Y = CY(IY)

            IF(IX .EQ. 1 ) THEN

              EY = CY(2)
              EZ = CY(3)
              IF(NRMA .EQ. 1) THEN
C--------------- Calculate t(Ey_max)...

                IF(EZ*EZ+EZ0*EZ0 .NE. 0.D0) THEN
C----------------- ...rather from t(Ey_max)  ==  t(Ez=0) ...
                  IF(EZ0*EZ .LE. 0.D0) THEN
                    TM = X - CX(1) / (EZ - EZ0) * EZ
                    JM = NOC
                  ENDIF

                ELSE
C----------------- ...or at worst from Ey[ t(Ey_max) ] ~ Ey_max
                  IF(EY .GT. 0.D0) THEN
                    IF(EY .LE. EY0) THEN
                      TM = X - CX(1)
                      JM = NOC
                    ENDIF
                  ENDIF
                ENDIF

              ELSEIF(NRMA .EQ. 2) THEN
              ENDIF

              IF(NRMC .GE. 1) THEN
C--------------- X-scale is omga_c * ( t - t(Ey_max) )
C                Y-scale is Ex,y,z / Ey_max

                IF(EY0*EY .LE. 0.D0) THEN
                  IF(EY .GT. EY0) THEN
                    TC1 = X - CX(1) / (EY - EY0) * EY
                    JC1 = NOC
                  ELSE
                    TC2 = X - CX(1) / (EY - EY0) * EY
                    JC2 = NOC
                  ENDIF
                  TC = .5D0 * (TC2 - TC1)                
                    if(tc.lt.0.D0) tc = -tc     !!!!! rustine
                ENDIF

                IF(CY(2) .GT. YM) YM = CY(2)

              ENDIF    ! NRMC

            ENDIF

            EY0 = EY
            EZ0 = EZ

          ELSEIF(NOC .EQ. 1) THEN

            X = 0.D0
            Y = CY(IY)

            EY0 = CY(2)
            EZ0 = CY(3)

          ENDIF        ! NOC

           CALL MINMAX(X,Y,
     >                     XMI,XMA,YMI,YMA)
           
          X = X + CX(IX)

        ELSEIF(KEP .EQ. 2) THEN
C----------- Plot

          IF(NOC .GE. 2) THEN

            Y = CY(IY)
            YDX = Y * CX(IX)
            XX = X

            IF(NRMA .EQ. 1) THEN
              IF(YM .NE. 0.D0) Y = Y / YM
              IF(IX .EQ. 1) THEN
C                XX = (X - TM) / TC 
                YDX = Y * CX(1) / TC
              ENDIF
            ELSEIF(NRMA .EQ. 2) THEN
            ENDIF

            IF(INRANG(XX,Y,XMI,XMA,YMI,YMA) ) THEN
              IF(NOC.EQ.NPTS) YDX = YDX/2.D0
              SYDX = SYDX + YDX
              CALL VECTPL(XX,Y,2) 
              IF(LIS .EQ. 2) CALL IMPV(NLOG,NOC,XX,Y,YDX,SYDX,IDUM)
            ELSE
              CALL VECTPL(XX,Y,4) 
            ENDIF

          ELSEIF(NOC .EQ. 1) THEN
            
            X = 0.D0
            Y = CY(IY)
            YDX = Y * CX(IX)

            IF(NRMA .EQ. 1) THEN
              IF(YM .NE. 0.D0) Y = Y / YM
              IF(IX .EQ. 1) THEN
                X = - TM / TC
                YDX = Y * CX(1) / TC
              ENDIF
            ELSEIF(NRMA .EQ. 2) THEN
            ENDIF

            CALL VECTPL(X,Y,4) 
            IF( INRANG(X,Y,XMI,XMA,YMI,YMA) ) THEN
C              SYDX = YDX 
              SYDX = YDX * 0.5D0

              CALL VECTPL(X,Y,2) 
              IF(LIS .EQ. 2) CALL IMPV(NLOG,NOC,X,Y,YDX,SYDX,IDUM)
            ELSE
              SYDX = 0.D0
            ENDIF 

          ENDIF   ! NOC = 1 

          X = X + CX(IX)

        ENDIF   ! KEP=1,2

      ELSEIF(NOC .EQ. 0) THEN
C------- Initialize scales

        XMI = 1.D10
        XMA = -1.D10
        YMI = 1.D10
        YMA = -1.D10
        IF(NRMC .GE. 1) YM = -1.D10

      ENDIF   ! NOC

      RETURN

      ENTRY SRPLTI
      NONEWT= .FALSE. 
C           call fbgtxt
C              write(*,*)  ' -----------NEWT-----------'
      RETURN
      END
