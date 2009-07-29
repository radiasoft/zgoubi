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
      SUBROUTINE PLTLPS(NLOG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH
      INCLUDE "MAXTRA.H"
      CHARACTER  LET
      COMMON/FAISC/ F(6,1),IMAX,LET(MXT),IEX(MXT),IREP(MXT)
      COMMON/OBJET/ FO(6,1),KOBJ,IDMAX,IMAXT
      COMMON/PTICUL/AM,Q,G,TOO
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DP   ,BR
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)
      COMMON/SYNCH/ RET(MXT), DPR(MXT)
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN
      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*7, KPOL(2)*9, KDIM(MXVAR)*7
      COMMON/INPVR/ KVAR, KPOL, KDIM
      COMMON/LUN/NDAT,NRES,NPLT,NFAI,NMAP,NSPN
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      LOGICAL OKPLOT
      SAVE OKPLOT
      DIMENSION X0(MXT),Y0(MXT),SX(MXT),SY(MXT)
      DIMENSION XM1(MXT), XM2(MXT),YM1(MXT),YM2(MXT)

C      DOUBLE PRECISION DUM

      IF(IPASS.EQ. 1 ) THEN
        OKPLOT = OKECH.AND.OKVAR
        IF(OKPLOT) THEN
          CALL LINTYP(-1)
        ELSE
          IF(.NOT.OKECH) WRITE(6,*)  ' PLEASE DEFINE SCALES '
          IF(.NOT.OKVAR) WRITE(6,*)  ' PLEASE DEFINE VARIABLES'
          STOP '  SBR PLTLPS'
        ENDIF
      ENDIF

C      CALL TXTFBG

      DO 1 I=1,IMAX
        IF(KX.LE. 5) THEN
          X=F(KX,I)
        ELSEIF(KX.EQ. 20) THEN
          X=IPASS
        ELSEIF(KX.EQ. 18) THEN
          X = RET(I)
        ELSEIF(KX.EQ. 19) THEN
          X=DPR(I) 
        ELSEIF(KX.LE. 24) THEN
          X=SF(KX-20,I)
        ELSEIF(KX.LE. 28) THEN
          X=SF(KX-24,I)
        ENDIF
        IF(KY.LE. 5) THEN
          Y=F(KY,I)
        ELSEIF(KY.EQ. 20) THEN
          Y=IPASS
        ELSEIF(KY.EQ. 18) THEN
          Y = RET(I)
        ELSEIF(KY.EQ. 19) THEN
          Y=DPR(I) 
        ELSEIF(KY.LE. 24) THEN
          Y=SF(KY-20,I)
        ELSEIF(KY.LE. 28) THEN
          Y=SF(KY-24,I)
        ENDIF

        IF(IPASS.EQ. 1)THEN
          X0(I)=X
          Y0(I)=Y
          SX(I) = 0.D0
          SY(I) = 0.D0
          XM1(I)=1.D10
          XM2(I)=-1.D10
          YM1(I)=1.D10
          YM2(I)=-1.D10
        ENDIF

        SX(I) = SX(I)+X
        SY(I) = SY(I)+Y
        IF(KY.GE. 25.AND.KY.LE. 28) Y=SY(I)/IPASS

        IF(OKPLOT) THEN
C          IF( IPASS.EQ. 1 .OR. (IPASS/(-KT))*(-KT) .EQ. IPASS) THEN
C           ... PLOT EACH KT TURN
            CALL VECTPL(X0(I),Y0(I),4)
            CALL VECTPL(X,Y,2)
            IF(LIS.EQ.2) 
     >           CALL IMPV(NLOG,IPASS,X,Y,DUM,DUM,IDUM)
C          ENDIF
        ENDIF

        CALL MINMAX(X,Y,XM1(I),XM2(I),YM1(I),YM2(I))
        X0(I)=X
        Y0(I)=Y

 1    CONTINUE       

C      CLOSE(NLOG)
      CALL FBGTXT

      IF(NRES.GT. 0) THEN
        WRITE(NRES,*)
        DO 2 I=1,IMAX
          WRITE(NRES,100) I,SX(I)/IPASS, SY(I)/IPASS
 100      FORMAT(20X,I3,' <X> = ',1PG12.4,' <Y> = ',1PG12.4)
          WRITE(NRES,101) XM1(I),XM2(I),YM1(I),YM2(I),X,Y
 101      FORMAT(23X,' XMI= ',1PG12.4,' XMA= ',1PG12.4
     >    ,/,23X,' YMI= ',1PG12.4,' YMA= ',1PG12.4
     >    ,/,' LAST   POINT :'
     >    ,/,23X,' X= ',1PG12.4,' Y= ',1PG12.4)
 2      CONTINUE       
      ENDIF
      
      RETURN
      END
