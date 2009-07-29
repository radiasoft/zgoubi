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
      SUBROUTINE LPSDRW(NLOG,YM,YPM,YMX,YPMX,U,A,B,KPS,XSIGU,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YM(*), YPM(*), U(*), A(*), B(*), XSIGU(*)
      DIMENSION YMX(*), YPMX(*)
      COMMON/CONST/ CL,PI,DPI,RAD,DEG,QE,AH
      LOGICAL OKECH, OKVAR, OKBIN 
      COMMON/ECHL/OKECH, OKVAR, OKBIN
      INCLUDE 'MXVAR.H'
      CHARACTER KVAR(MXVAR)*7, KPOL(2)*9, KDIM(MXVAR)*7
      COMMON/INPVR/ KVAR, KPOL, KDIM
      INCLUDE 'MAXNTR.H'             
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB 

      CHARACTER REP
        
      CHARACTER HVL(3)*12, TXT*80
      DATA HVL / 'Horizontal', 'Vertical', 'Longitudinal' /

      SAVE JJ
      SAVE BET,ALP,EPS,CLIP,CLIPP

      LOGICAL AGAIN
      DATA AGAIN / .FALSE. /

      GOTO 21

 20   CONTINUE      
      CALL FBGTXT
      WRITE(*,FMT='(/,''  Press RETURN for more'')') 
      READ(*,200,ERR=20) REP 
 200  FORMAT(A1)

 21   CONTINUE
      CALL FBGTXT
      CALL HOMCLR

      WRITE(*,104) 
 104  FORMAT(5X,' OPTIONS : ')
      IF(KPS.EQ. 0) WRITE(*,109) '(Yo, Zo, Do)' 
      IF(KPS.EQ. 1) WRITE(*,109) '(Y, Z, D)' 
 109  FORMAT           
     1(/,5X,' 1, 2, 3  ** Plot the fitting ellipse ',A,' **')
      WRITE(*,110) XSIGU(1),XSIGU(2),XSIGU(3)
 110  FORMAT(
     4   5X,'    4  ** Plot the phase-space points sample **'
     5,/,5X,'    5  H/V borders at ',F5.2,'/',
     >                      F5.2,'/',F5.2,' * sigma(Emittance)'
     6,/,5X,'    6  Plot an ellipse '
     8,/,5X,'    8  Print screen'
     9,/,5X,'    9  EXIT  THIS  MENU '
     2,/,5X,'   12  ERASE  DISPLAY    '
     >,/)

      WRITE(*,100)
 100  FORMAT('$  Option  number : ',$)
      READ(*,101,ERR=21) IOPT
101   FORMAT(I2)
      GOTO ( 1, 1, 1, 4, 5, 6,21, 8, 9,21,21,12) IOPT  
      GOTO 21


 1    CONTINUE
      KX0 = KX
      KY0 = KY
      IF(IOPT.EQ. 1) THEN
        JJ=1
        KX=2
        IF(KPS.EQ. 0) KX=12
        KY=KX+1
      ELSEIF(IOPT.EQ. 2) THEN
        JJ=2
        KX=4
        IF(KPS.EQ. 0) KX=14
        KY=KX+1
      ELSEIF(IOPT.EQ. 3) THEN
        JJ=3
        KX=7   !6  !18
        KY=20  !1  !19
        IF(KPS.EQ. 0) THEN
          KX=16
          KY=11
        ENDIF
      ENDIF

      OKVAR = .TRUE.

      IF(U(JJ).EQ. 0.D0) THEN
        WRITE(*,131) HVL(JJ)
 131    FORMAT(/,10X,A12,'  PHASE-SPACE  ELLIPSE'
     >  ,/,'  EMITTANCE = 0.'
     >  ,/,'  ELLIPSE  PARAMETERS  UNDETERMINED'
     >  ,/,'  NO  PLOT  POSSIBLE'
     >  ,/,'  DID  YOU  FIT,  FIRST ? - IF NOT, GO TO PREVIOUS MENU',/)

      ELSE      
C        CALL TXTFBG
        IF(.NOT. OKECH) THEN 
          J = 2*JJ - 1
          XMI = YMX(J)
          XMA = YMX(J+1)
          YMI = YPMX(J)
          YMA = YPMX(J+1)
          IF(XMI.LT. XMA .AND. YMI.LT. YMA) THEN
C            CALL TXTFBG        
            OKECH=.TRUE.
          ELSE
            WRITE(*,140) HVL(JJ)
 140        FORMAT(/,10X,A12,'  PLOT OF PHASE-SPACE POINTS NOT POSSIBLE'
     >        ,/,'  Min-max problem (Presumably zero emittance)',/)
          ENDIF
        ENDIF

        IF(OKECH) THEN
          CALL TRAXES(XMI,XMA,YMI,YMA,1)
          CALL LINTYP(1)
C          R=XSIGU(JJ) * SQRT(B(JJ)*U(JJ))
          R= SQRT(XSIGU(JJ) * U(JJ) * B(JJ)) 
          Y = R + YM(JJ)
C FM 03/02
C          YP= A(JJ)/B(JJ) * R + YPM(JJ)
          YP= -A(JJ)/B(JJ) * R + YPM(JJ)
CC------------ MAY 1999
C               Y = Y /2.D0
C               YP = YP /2.D0
          CALL VECTPL(Y,YP,4)      
          DO 201 IPHI=0,360
            PHI = IPHI*DPI/360
            Y = R * COS(PHI)
C            YP= -(A(JJ)* Y + R * SIN(PHI) )/B(JJ)  + YPM(JJ)
C FM 03/02
C            YP= (A(JJ)* Y + R * SIN(PHI) )/B(JJ)  + YPM(JJ)
            YP= (-A(JJ)* Y + R * SIN(PHI) )/B(JJ)  + YPM(JJ)
            Y = Y + YM(JJ)
CC------------ MAY 1999
C               Y = Y /2.D0
C               YP = YP /2.D0
            CALL VECTPL(Y,YP,2)      
            IF(LIS .EQ. 2)  CALL IMPV(NLOG,0,Y,YP,DUM,DUM,IDUM)
 201      CONTINUE

          CALL FBGTXT
          WRITE(*,108) KVAR(KY),KVAR(KX)
 108      FORMAT(A,' v.s. ',A) 
          J = 2*JJ - 1
          WRITE(*,107) YMX(J),YMX(J+1),YPMX(J),YPMX(J+1)
 107      FORMAT(' Min-max - Hor.: ',1P,2G13.5,'; Ver.: ',2G13.5)
          WRITE(*,105) U(JJ),B(JJ),A(JJ),YM(JJ),YPM(JJ)
 105      FORMAT(' Eps/pi, Beta, Alpha, center : ',1P,
     >                                G12.4,4G11.3,' (MKSA)')

          WRITE(TXT,108) KVAR(KY),KVAR(KX)
          CALL TRTXT(10D0,26.D0,TXT,80,0)
          WRITE(TXT,107) YMX(J),YMX(J+1),YPMX(J),YPMX(J+1)
          CALL TRTXT(10.D0,10.D0,TXT,80,0)
          WRITE(TXT,105) U(JJ),B(JJ),A(JJ)
          CALL TRTXT(10.D0,18.D0,TXT,80,0)

          CALL TRTXT(0.0D0,0.0D0,' ',45,0)

          CALL LINTYP(-1)

        ENDIF

      ENDIF

      GOTO 20

 4    CONTINUE
        J=2*JJ-1
        IF(.NOT. OKECH) THEN 
          J = 2*JJ - 1
          XMI = YMX(J)
          XMA = YMX(J+1)
          YMI = YPMX(J)
          YMA = YPMX(J+1)
          IF(XMI.LT. XMA .AND. YMI.LT. YMA) THEN
            CALL TRAXES(XMI,XMA,YMI,YMA,1)
            OKECH=.TRUE.
          ELSE
            WRITE(*,140) HVL(JJ)
          ENDIF
        ENDIF
        IF(OKECH) THEN
          CALL LINTYP(9)
          CALL VECTPL(COOR(1,J),COOR(1,J+1),4)
          NPOINT=NPTS
C          IF(NPOINT.GT.5000) THEN
          IF(NPOINT.GT.NTRMAX) THEN
C            NPOINT=5000
            NPOINT=NTRMAX
            WRITE(6,*) 
            WRITE(6,*) ' *** Warning from subroutine lpsdrw : '
            WRITE(6,*) '     Storage capacity limited to ',NTRMAX,
     >                  ' particles'
          ENDIF
          DO 42 I=1,NPOINT
            CALL VECTPL(COOR(I,J),COOR(I,J+1),2)
            IF(LIS.EQ.2) 
     >        CALL IMPV(NLOG,I,COOR(I,J),COOR(I,J+1),DUM,DUM,IDUM)
 42       CONTINUE     
          CALL LINTYP(-1)
          CALL FBGTXT
          CALL LPSCNT(YM,YPM,U,A,B,XSIGU,NLOG,
     >                                        NCOUNT)
        ENDIF
      GOTO 20

 5    CONTINUE
        WRITE(*,*) ' Give ellipse extents in units of sigma'
        WRITE(*,*) '   [Now :',XSIGU(1), XSIGU(2), XSIGU(3),']  :'
        READ(*,*,ERR=5) XSIGU(1),  XSIGU(2),  XSIGU(3)
      GOTO 21

 6    CONTINUE
        WRITE(*,FMT='(''Give ellipse parameters and center'')')
        WRITE(*,61) EPS,BET,ALP,CLIP,CLIPP
 61     FORMAT('  (Now  EPS/pi BET ALP x xp = ',1P,5G12.4,')   : ',$)
C        READ(*,*,ERR=6) EPS,BET,ALP,CLIP,CLIPP
        READ(*,*,ERR=62) EPS,BET,ALP,CLIP,CLIPP
        AGAIN = .FALSE.
        GOTO 63
 62     CONTINUE
        AGAIN = .TRUE.
        READ(88,*,ERR=6,END=6) CMAX, EPS,ALP,BET,CLIP,CLIPP
 63     CONTINUE
        IF(.NOT.OKECH) THEN
          R= SQRT(EPS*BET)
          Y = R + CLIP
          YP= -ALP/BET * R + CLIPP
          XMI = 1.D10
          XMA = -1.D10
          YMI = 1.D10
          YMA = -1.D10
          DO 602 IPHI=0,360
            PHI = IPHI*DPI/360
            Y = R * COS(PHI)
            YP= (-ALP* Y + R * SIN(PHI) )/BET  + CLIPP
            Y = Y + CLIP
            IF(Y.LT.XMI) XMI = Y
            IF(Y.GT.XMA) XMA = Y
            IF(YP.LT.YMI) YMI = YP
            IF(YP.GT.YMA) YMA = YP
 602      CONTINUE
          IF(XMI.LT. XMA .AND. YMI.LT. YMA) THEN
            OKECH=.TRUE.
            CALL TRAXES(XMI,XMA,YMI,YMA,1)
          ELSE
            WRITE(*,140) HVL(JJ)
          ENDIF
        ENDIF
        CALL LINTYP(1)
        R= SQRT(EPS*BET)
        Y = R + CLIP
        YP= -ALP/BET * R + CLIPP
        CALL VECTPL(Y,YP,4)      
        DO 601 IPHI=0,360
          PHI = IPHI*DPI/360
          Y = R * COS(PHI)
          YP= (-ALP* Y + R * SIN(PHI) )/BET  + CLIPP
          Y = Y + CLIP
          CALL VECTPL(Y,YP,2)      
          IF(LIS .EQ. 2)  CALL IMPV(NLOG,0,Y,YP,DUM,DUM,IDUM)
 601    CONTINUE
        CALL LINTYP(-1)
        CALL FBGTXT
        YM(JJ)=CLIP
        YPM(JJ)=CLIPP
C        XSIGU(JJ)=EPS/U(JJ)
        U(JJ)=EPS
        A(JJ)=ALP
        B(JJ)=BET
        CALL LPSCNT(YM,YPM,U,A,B,XSIGU,NLOG,
     >                                      NCOUNT)
      IF(AGAIN) GOTO 6
      GOTO 20

 8    CONTINUE
C        CALL MENVCF
        CALL SAVPLT
      GOTO 21

 9    CONTINUE
      CALL LINTYP(-1)
      RETURN 1

 12   CONTINUE
      CALL CLSCR
      OKECH=.FALSE.
      GOTO 21

      END
