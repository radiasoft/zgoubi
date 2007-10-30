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
      SUBROUTINE PLHIST(NL,LM,NPTS,NPTR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXB=1000)
      COMMON/B/BINS(MXB+2)
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      LOGICAL OKBIN,OKECH,OKXAV,OKYAV

      OKECH = .TRUE.
      CALL BIN3W(0,0)

C------ Histogram X
      CALL BIN(NL,LM,OKECH,KX,NB,LIS,NPTS,NPTR,
     >                        NOCE,OKBIN,XMOY,SIG,XMI,XMA,*97)
        
        BMI = 1.D10
        BMA =-1.D10
        DO 111 I=1,NB
          Y=BINS(I)
          IF(Y .GT. BMA) BMA=Y
          IF(Y .LT. BMI) BMI=Y
 111    CONTINUE
        XMIB=BINS(MXB+1)
        XMAB=XMIB + (NB-1.D0)*BINS(MXB+2)

        CALL PLOT4R(OKXAV,OKYAV)
        IF(OKXAV) THEN
          DO I=1,NB
            BINS(NB) = BINS(NB) - XMOY
          ENDDO
          BINS(MXB+1) = BINS(MXB+1) - XMOY
        ENDIF

      CALL LINTYP(1)
C A l'interieur de la fenetre
C      FB = -(YMA - YMI)/10.D0
      FB = -(YMA - YMI)/4.D0
C A l'exterieur de la fenetre... 
C      FB = (YMA - YMI)/10.D0
      SYDX=0.D0
      YMAX = -1.D10
      DB=BINS(MXB+2)
      X=BINS(MXB+1)
      CALL VECTPL(X,YMA,4)
      DO 1 IB=1,NB
        Y=BINS(IB) 
        YY = Y/(BMA-BMI) * FB + YMA
        CALL VECTPL(X,YY,2)
        SYDX = SYDX + Y*DB
        IF(LIS .EQ. 2) CALL IMPV(NLOG,NB,X,YY,Y*DB,SYDX,IDUM,IDUM,KX,KY)
        X=X+DB
        CALL VECTPL(X,YY,2)
        IF(Y.GT.YMAX) THEN
          YMAX=Y
          XYMAX=X
        ENDIF
 1    CONTINUE
      CALL VECTPL(X,YMA,2)
      CALL FBGTXT

 97   CONTINUE

C------ Histogram Y
      CALL BIN(NL,LM,OKECH,KY,NB,LIS,NPTS,NPTR,
     >                        NOCE,OKBIN,XMOY,SIG,YMI,YMA,*98)
        BMI = 1.D10
        BMA =-1.D10
        DO 112 I=1,NB
          Y=BINS(I)
          IF(Y .GT. BMA) BMA=Y
          IF(Y .LT. BMI) BMI=Y
 112   CONTINUE
        XMIB=BINS(MXB+1)
        XMAB=XMIB + (NB-1.D0)*BINS(MXB+2)

        IF(OKYAV) THEN
          DO I=1,NB
            BINS(NB) = BINS(NB) - XMOY
          ENDDO
          BINS(MXB+1) = BINS(MXB+1) - XMOY
        ENDIF

      CALL LINTYP(1)
C      FB = -(XMA - XMI)/10.D0
      FB = -(XMA - XMI)/4.D0
C      FB = (XMA - XMI)/10.D0
      SYDX=0.D0
      YMAX = -1.D10
      DB=BINS(MXB+2)
      Y=BINS(MXB+1)
      CALL VECTPL(XMA,Y,4)
      DO 2 IB=1,NB
        X=BINS(IB) 
        XX = X/(BMA-BMI) * FB + XMA
        CALL VECTPL(XX,Y,2)
        SYDX = SYDX + Y*DB
        IF(LIS .EQ. 2) CALL IMPV(NLOG,NB,XX,Y,Y*DB,SYDX,IDUM,IDUM,KX,KY)
        Y=Y+DB
        CALL VECTPL(XX,Y,2)
        IF(Y.GT.YMAX) THEN
          YMAX=X
          XYMAX=Y
        ENDIF
 2    CONTINUE
      CALL VECTPL(XMA,Y,2)
      CALL FBGTXT
      CALL LINTYP(-1)
      RETURN
 98   WRITE(6,*)' *** In SBR PLHIST : Error Xmi-max while binning histo'
      RETURN
      END
