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
      SUBROUTINE PLTHIS(NLOG,NL,NOC,OKBIN,OKECH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKBIN, OKECH
C     ----------------
C     PLOT HISTOGRAMME
C     ----------------
      PARAMETER (MXB=1000)
      COMMON/B/BINS(MXB+2)
      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*7, KPOL(2)*9, KDIM(MXVAR)*7
      COMMON/INPVR/ KVAR, KPOL, KDIM
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      CHARACTER TXT*80

      IF(OKBIN) THEN
        CALL BIN2(NOC,
     >               XMOY,SIG)
      ELSE
        CALL FBGTXT
        CALL BIN(NL,OKECH,KX,NB,
     >                          NOC,OKBIN,XMOY,SIG,XMI,XMA,*99)
      ENDIF

      CALL LINTYP(1)

      SYDX=0.D0
      YMAX = -1.D10
      DB=BINS(MXB+2)
      X=BINS(MXB+1)
      CALL VECTPL(X,0.D0,4)
      DO 1 IB=1,NB
        Y=BINS(IB)
        CALL VECTPL(X,Y,2)
        SYDX = SYDX + Y*DB
        IF(LIS .EQ. 2) CALL IMPV(NLOG,NB,X,Y,Y*DB,SYDX,IDUM)
C        if(ib.gt.1) write(78,*) x,-BINS(IB)+BINS(IB-1), '   sbr plthis '
        X=X+DB
        CALL VECTPL(X,Y,2)
        IF(Y.GT.YMAX) THEN
          YMAX=Y
          XYMAX=X
        ENDIF
 1    CONTINUE
      CALL VECTPL(X,0.D0,2)
      CALL LINTYP(-1)

      CALL TRKVAR(NB,KVAR(KY),KDIM(KY),KVAR(KX),KDIM(KX)) 

      WRITE(6,100)
 100  FORMAT('Mean',12X,'Sigma',12X,'X(max)',12X,'Counts')
      WRITE(6,101) XMOY,SIG,XYMAX,NOC
 101  FORMAT(1P,3G16.9,I7)
      WRITE(TXT,102) XMOY,SIG,XYMAX,NOC
 102  FORMAT('Avrg:',1P,G11.5,'; Sgma:',G11.5,'; X(max):',G11.5,
     > '; Counts:',I7)
      CALL TRTXT(10.D0,21.D0,TXT,0)

      RETURN
 99   WRITE(6,*) ' SBR PLTHIS   ERROR DURING BINING'
      END
