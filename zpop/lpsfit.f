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
      SUBROUTINE LPSFIT(NLOG,KPR,LM,
     >                              YM,YPM,YMX,YPMX,U,A,B,*,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YM(*), YPM(*), U(*), A(*), B(*)
      DIMENSION YMX(*), YPMX(*)
      INCLUDE 'MAXNTR.H'
      COMMON/TRACKM/COOR(NTRMAX,9),NPTS,NPTR

      DIMENSION G(3)
      CHARACTER TXT(3)*5, REP*1

      INCLUDE 'MAXCOO.H'

      DIMENSION UM(3), UMI(3), UMA(3), SMEAR(3)
      CHARACTER*1 KLET

      PARAMETER ( PI=3.1415926536, SQ2 = 1.414213562)
C      PARAMETER ( PI=4.D0*ATAN(1.D0), SQ2 = SQRT(2.D0) )

      DATA TXT/ 'HORI.', 'VERT.', 'LONG.'/

      DO 7 J=1,MXJ-1,2
        YMX(J) = 1.D10
        YMX(J+1) = -1.D10
        YPMX(J) = 1.D10
        YPMX(J+1) = -1.D10
 7    CONTINUE

      DO 2 J=1,MXJ-1,2
C------- JJ = 1 , 2 or 3  for  Y-T, Z-P or T-P(time-momentum) planes
        J1 = J
        J2 = J+1
        JJ= J2 / 2
        XM=0.D0
        XPM=0.D0
        SNPT = 0.D0
        DO 21 I=1,NPTS
            X = COOR(I,J1)
            XP = COOR(I,J2)
            XM = XM + X
            XPM = XPM + XP
c              write(*,fmt='(a,2g12.4,3i4)') 'blibli ',x,xp,i,npts,jj
C----------- Min-max of the distribution
            IF(YMX(J1) .GT. X) YMX(J1) = X
            IF(YMX(J2) .LT. X) YMX(J2) = X
            IF(YPMX(J1) .GT. XP) YPMX(J1) = XP
            IF(YPMX(J2) .LT. XP) YPMX(J2) = XP
            SNPT = SNPT + 1
 21     CONTINUE
        XM = XM/SNPT
        XPM = XPM/SNPT
        YM(JJ)=XM
        YPM(JJ)=XPM
 2    CONTINUE

      DO 25 J=1,MXJ-1,2
C------- JJ = 1 , 2 or 3  for  Y-T, Z-P or T-P(time-momentum) planes
        J1 = J
        J2 = J+1
        JJ= (J+1) / 2
        X2=0.D0
        XP2=0.D0
        XXP=0.D0
        SNPT = 0.D0
        DO 26 I=1,NPTS
            X = COOR(I,J1)
            XP = COOR(I,J2)
            X2 = X2 + (X-YM(JJ))**2
            XP2 = XP2 + (XP-YPM(JJ))**2
            XXP = XXP + (X-YM(JJ))*(XP-YPM(JJ))
            SNPT = SNPT + 1
c              write(*,fmt='(a,3g12.4,i4)') 'blublu ',x,x2,ym(jj),jj
 26     CONTINUE
        X2  = X2/SNPT
        XP2 = XP2/SNPT
        XXP = XXP/SNPT

C G. Leleux : surface de l'ellipse S=4.pi.sqrt(DELTA)
C Soit d11=X2/sqrt(DELTA), d12=XXP/sqrt(DELTA), d22=XP2/sqrt(DELTA), alors 
C  d22.x^2-2.d12.x.x'+d11.x'^2=S/pi=4sqrt(DELTA), ce qui permet d'ecrire 
C   gamma=d22=XP2/sqrt(DELTA), -alpha=d12=XXP/sqrt(DELTA), beta=d11=X2/sqrt(DELTA). 
C En outre, par definition des dij, 
C     2.sigma_x=sqrt(d11.S/pi),  2.sigma_x'=sqrt(d22.S/pi). 
C En outre, frontiere : 
C          <x^2>_frontiere=2.(sigma_x)^2,    <x'^2>_frontiere=2.(sigma_x')^2

        SQ = SQRT(X2*XP2-XXP*XXP) 
        IF(SQ .GT. 0.D0) THEN
          B(JJ)=  X2/SQ
C Error  FM 03/02
C          A(JJ)=  XXP/SQ
          A(JJ)=  -XXP/SQ
          G(JJ)=  XP2/SQ
        ENDIF
C------- Courant invariant at 1 sigma is U=4.sqrt(DELTA)=Eps/pi (consistant with zgoubi !!) :
C Eps=ellipse surface
C        U(JJ) = 4.D0*SQ
        U(JJ) = SQ

c        write(*,fmt='(a,3g12.4,2i6)') 'blabla ',x2,xp2,u(jj),npts,jj

 25   CONTINUE

      IF(KPR .EQ. 0) RETURN
     
C----- SMEAR
      DO 3 J=1,MXJ-1,2
C------- JJ = 1 , 2 or 3  for  Y-T, Z-P or T-P(time-momentum) planes
        J1 = J
        J2 = J+1
        JJ= (J+1) / 2
        UMA(JJ) = -1.D10
        UMI(JJ) = 1.D10
        UM(JJ)=0.D0
        U2M=0.D0
        SNPT = 0.D0
        DO 31 I=1,NPTS
C--------- Normalized coordinates (*Beta), for phase-space point I:
            X = ( COOR(I,J1) - YM(JJ) )
            XP = ( COOR(I,J2) - YPM(JJ) )
            XN = X 
            XPN = ( A(JJ) * XN + B(JJ) * XP )
C----------- Courant invariant Epsilon/pi at phase-space point I:
            UI = ( XN * XN + XPN * XPN )/ B(JJ)

            IF(UI .GT. UMA(JJ)) UMA(JJ) = UI
            IF(UI .LT. UMI(JJ)) UMI(JJ) = UI

            UM(JJ) = UM(JJ) + UI
            U2M = U2M + UI * UI
            SNPT = SNPT + 1
 31     CONTINUE

        UM(JJ) = UM(JJ)/SNPT
        U2M = U2M/SNPT
        SMEAR(JJ) = SQRT( U2M - UM(JJ) * UM(JJ) )

 3    CONTINUE

C----- Twiss parameters and emittance
      CALL READC5(NT1,NT2)
      CALL READC9(KKEX,KLET)
      WRITE(*,FMT='(/,'' Particle # '',I6,'' to '',I6,''    @ lmnt # '',
     >I5,''. '',I6,'' points,  IEX='',I2)') NT1,NT2, LM, NINT(SNPT),KKEX
      WRITE(*,100) 
 100  FORMAT(/,10X,'Ellipse parameters and center  (MKSA) :',//,T12,
     >'BETA',T24,'ALPHA',T36,'GAMMA',T48,'EPS/PI',T63,'CENTER')

      DO 6 JJ=1,3
        I=2*JJ-1
C        IF(U(JJ).LE. 1.D-30) THEN
        IF(U(JJ).EQ. 0.D0) THEN
          WRITE(*,121) TXT(JJ),U(JJ)
 121      FORMAT(A,T20,'*** UNDETERMINED ***',T49,1P,G12.4)
        ELSE
          WRITE(*,120) TXT(JJ),B(JJ),A(JJ),G(JJ),U(JJ),YM(JJ),YPM(JJ)
 120      FORMAT(A,T10,1P,G12.4,T22,G12.4,T34,G12.4, 
     >    T46,G12.4,T58,2G12.4)
        ENDIF
 6    CONTINUE

      WRITE(*,FMT='(/,A)') ' Frontier values : '
      DO 66 JJ=1,3
        I=2*JJ-1
        FAC = 1.D0
        IF(U(JJ).EQ. 0.D0) THEN
          WRITE(*,121) TXT(JJ),U(JJ)
        ELSE
          WRITE(*,120) TXT(JJ),B(JJ),A(JJ),G(JJ),U(JJ),YM(JJ),YPM(JJ)
        ENDIF
 66   CONTINUE

C----- Smear
      WRITE(*,130) 
 130  FORMAT(/,10X,'Smear :',/,T17,
     >'<E/pi>',T26,'Sigma(E/pi)/<E/pi>',T48,'Min(E/pi)',T60,'Max(E/pi)')

      DO 5 JJ=1,3
        I=2*JJ-1
C        IF(U(JJ).LE. 1.D-30) THEN
        IF(U(JJ).EQ. 0.D0) THEN
          WRITE(*,131) TXT(JJ),U(JJ)
 131      FORMAT(A,T13,1P,G12.4,'      *** UNDETERMINED ***')
        ELSE
          WRITE(*,132) TXT(JJ),U(JJ),SMEAR(JJ)/U(JJ),UMI(JJ),UMA(JJ)
 132      FORMAT(A,T13,1P,G12.4,T28,G12.4,T45,G12.4,T57,G12.4)
        ENDIF
 5    CONTINUE

      IF(KPR.EQ.1) THEN   
 20     WRITE(*,*)
        WRITE(*,*) '  PRINT IN zpop.log (Y/N)?' 
        READ(*,FMT='(A1)',ERR=20) REP        
        IF(REP .NE. 'N' .AND. REP .NE. 'n') REP = 'y'
      ELSEIF(KPR.EQ.2) THEN
        REP = 'y'
      ENDIF

      IF(REP.EQ. 'y') THEN

C----- Twiss parameters and emittance
        WRITE(*,*) '  Ellipse parameter calculated will'
     >  ,' be printed in zpop.log'

        CALL READC9(KKEX,KLET)
        WRITE(NLOG,FMT='(/,'' Particle # '',I6,'' to '',I6,
     >  ''  @ lmnt # '',I5,''. '',I6,'' points,  IEX='',I2)') 
     >  NT1,NT2,LM,NINT(SNPT),KKEX
        WRITE(NLOG,100) 

        DO 4 JJ=1,3
          IF(U(JJ).LE. 1.D-15) THEN
            WRITE(NLOG,121) TXT(JJ),U(JJ) 
C           WRITE(*,121) TXT(JJ),U(JJ) 
          ELSE
            WRITE(NLOG,120) 
     >         TXT(JJ),B(JJ),A(JJ),G(JJ),U(JJ),YM(JJ),YPM(JJ)
          ENDIF
 4      CONTINUE

C------- Smear
        WRITE(NLOG,130) 
        DO 51 JJ=1,3
          I=2*JJ-1
          WRITE(NLOG,131) TXT(JJ),U(JJ),SMEAR(JJ),UMI(JJ),UMA(JJ)
 51     CONTINUE

        RETURN 1
      ENDIF      

      RETURN 2
      END
