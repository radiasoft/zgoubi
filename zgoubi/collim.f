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
C  USA
C  -------
      SUBROUTINE COLLIM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -------------------------------------------------------
C     COLLIM EST EQUIVALENT A CHAMBR DE DIMENSION LONGITUDI-
C     NALE NULLE , ET A LES MEMES EFFETS ( UPDATES SORT(I) ET
C     FAIT IEX(I)=-4 SI UNE PARTICULE SORT )
C     -------------------------------------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH
 
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER(1) LET
      INCLUDE "C.FAISCT.H"     ! COMMON/FAISCT/ LET(MXT)
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)
  
C      DIMENSION YM(MXJ), YPM(MXJ), U(MXJ), ALPHA(MXJ), BETA(MXJ)
C      DIMENSION YMX(MXJ), YPMX(MXJ)

      IL =   NINT(A(NOEL,1))
      IFRM = NINT(A(NOEL,10))
      JFRM = NINT(10.D0*A(NOEL,10)) - 10*IFRM
      A1 = A(NOEL,11)
      A2 = A(NOEL,12)
      A3 = A(NOEL,13)
      A4 = A(NOEL,14)
      XM = A(NOEL,15)
      XPM = A(NOEL,16)

      IF(IL .EQ. 0) THEN

        IF(NRES .GT. 0) WRITE(NRES,100)
 100    FORMAT(/,20X,' +++++++  Collimator  is  inactive  +++++++',/)

        RETURN
      ENDIF

      IF(IFRM .LE. 2) THEN
        IF(JFRM.EQ.1) THEN
          YL=.5D0*(A2-A1)
          ZL=.5D0*(A4-A3)
          YC=.5D0*(A1+A2)
          ZC=.5D0*(A3+A4)
        ELSE
          YL=A1
          ZL=A2
          YC=A3
          ZC=A4
        ENDIF
      ELSE
        YL = A1
        ZL = A2
        YC = A3
        ZC = A4
        IF(IFRM .GE. 11) THEN
C--------- Phase-space (acceptance) collimator
          EPSPI = YC
          CUTOFF = ZC

          IF(IFRM .LE. 13) THEN
            ALP = YL
            BET = ZL
          ELSE
              IF(IFRM .EQ. 14) THEN
C--------------- Compute H matched llips
                JJ = 1
                CALL LPSFIT(JJ, 
     >                                   EMIT,ALP,BET,XM,XPM)
              ELSEIF(IFRM .EQ. 15) THEN
C--------------- Compute V matched llips
                JJ = 2
                CALL LPSFIT(JJ, 
     >                                   EMIT,ALP,BET,XM,XPM)
              ELSEIF(IFRM .EQ. 16) THEN
C---------------Compute LONGITUDINAL matched llips
                JJ = 3
                CALL LPSFIT(JJ, 
     >                                   EMIT,ALP,BET,XM,XPM)
              ENDIF
          ENDIF  
        ENDIF
      ENDIF

        IF(NRES .GT. 0) THEN
          IF    (IFRM .EQ. 1) THEN
            WRITE(NRES,101)
 101        FORMAT(/,20X,' Rectangular collimator.  ')
            WRITE(NRES,103) YL,ZL,YC,ZC
 103        FORMAT(25X,' YL = +/-',1P,G12.4,' cm'
     >             ,5X,' ZL = +/-',   G12.4,' cm'
     >       ,//,20X,' Centred  at :',/
     >            ,25X,' YC =    ',   G12.4,' cm'
     >             ,5X,' ZC =    ',   G12.4,' cm')
            WRITE(NRES,107) YC-YL,YC+YL,ZC-ZL,ZC+ZL
 107        FORMAT(/,20X,' => Max coordinates :    ',1P,G12.4,' < Y < '
     >      ,G12.4,'  CM',/,45X,G12.4,' < Z < ',G12.4,'  CM')
          ELSEIF(IFRM .EQ. 2) THEN
            WRITE(NRES,102)
 102        FORMAT(/,20X,' Elliptical collimator.  ')
            WRITE(NRES,103) YL,ZL,YC,ZC
            WRITE(NRES,108) YC,YL,ZC,ZL
 108        FORMAT(/,20X,' => Max coordinates :',/,25X,1P,'((Y-',G10.3
     >      ,')/',G10.3,')**2 + ((Z-',G10.3,')/',G10.3,')**2 < 1')
          ELSEIF(IFRM .GE. 11) THEN
C------------- Phase-space (acceptance) collimator
            IF(IFRM .EQ.11 .OR. IFRM .EQ.14.OR. IFRM .EQ.17) THEN
              
              WRITE(NRES,FMT='(/,20X,
     >            '' Horizontal  acceptance  collimator.'')')
            ELSEIF(IFRM .EQ.12 .OR. IFRM .EQ.15.OR. IFRM .EQ.18) THEN
              WRITE(NRES,FMT='(/,20X,
     >            '' Vertical  acceptance  collimator.'')')
            ELSEIF(IFRM .EQ.13 .OR. IFRM .EQ.16.OR. IFRM .EQ.19) THEN
              WRITE(NRES,FMT='(/,20X,
     >            '' Longitudinal  acceptance  collimator.'')')
            ENDIF
            WRITE(NRES,FMT='(/,15X,''     Beta,      Alpha,'',
     >      ''    E/pi,     Cut-off (units of Bet*Eps/pi),  XM, XPM :'',
     >      /,T16,1P,3G12.4,A6,G10.2,A2,15X,2G12.4)') 
     >      BET,ALP,EPSPI,'     (',CUTOFF,' )',XM, XPM
          ENDIF
        ENDIF

      LIVE = 0
      N1 = 0
      DO 1 I=1,IMAX
        IF(IEX(I) .LT. -1) GOTO 1
        LIVE = LIVE + 1
  
        IF    (IFRM .LT. 10) THEN

          IF    (IFRM .LT. 6) THEN
C----------- Physical collimator

            YP2 = (F(2,I) - YC)/YL
            YP2 = YP2*YP2
            ZP2 = (F(4,I) - ZC)/ZL
            ZP2 = ZP2*ZP2
 
            IF    (IFRM .EQ. 1) THEN
C-------------- Rectangular
 
              IF( YP2 .GE. 1.D0 .OR. ZP2 .GE. 1.D0 ) THEN
                N1 = N1 + 1
C---------------- SORT = path length of a particle when stopped 
                SORT(I) = F(6,I)
                IF(NRES .GT. 0) THEN
                  IF(IL .EQ. 2) THEN
                    IF(N1 .EQ. 1) WRITE(NRES,104)
 104                FORMAT(5X,/,'Coordinates  of  particles  stopped :'
     >              ,//,7X,'NTOT',T40,'(Do,Yo,To,Zo,Po)',
     >              T93,'(D,Y,T,Z,P)',/)
                    CALL CNTOUR(
     >                          NOUT)
                    WRITE(NRES,106) N1,NOUT,LET(I),IEX(I), 
     >               (FO(J,I),J=1,5),(F(J,I),J=1,5),I,IREP(I)
 106                FORMAT(1X,I3,I7,2X,A1,2X,I2,2(F8.4,4F10.3,2X),2I3)
                    WRITE(NRES,105) F(6,I)
 105                FORMAT(10X,'  path  length  =',F10.2,' cm')
                  ENDIF
                ENDIF
                CALL KSTOP(4,I,IEX(I),*10)
 10             CONTINUE
              ENDIF
 
            ELSEIF(IFRM .EQ. 2) THEN
C-------------- Elliptical
 
              IF( ( YP2 + ZP2 ) .GE. 1.D0 ) THEN
                N1 = N1 + 1
                SORT(I) = F(6,I)
                IF(NRES .GT. 0) THEN
                  IF(IL .EQ. 2) THEN
                    IF(N1 .EQ. 1) WRITE(NRES,104)
                    CALL CNTOUR(
     >                          NOUT)
                    WRITE(NRES,106) N1,NOUT,LET(I),IEX(I),
     >               (FO(J,I),J=1,5),(F(J,I),J=1,5),I,IREP(I)
                    WRITE(NRES,105) F(6,I)
                  ENDIF
                ENDIF
                CALL KSTOP(4,I,IEX(I),*11)
 11             CONTINUE
              ENDIF
 
            ENDIF

          ELSEIF(IFRM .GE. 6) THEN
C------------ Longitudinal space  collimator
C            IFRM = 6, 7 for rspctvly S, Time
C            JFRM = 1, 2 for rspctvly 1+dp/p, kinE
            JH = IFRM
            JV = JFRM
            H = F(JH,I)
            IF    (JFRM.EQ.1) THEN
              V = F(JV,I)
            ELSEIF(JFRM.EQ.2) THEN
              P = BORO*CL9 *F(1,I) * AMQ(2,I)
              ENERG = SQRT(P*P + AMQ(1,I)*AMQ(1,I))
              V = ENERG - AMQ(1,I)
            ENDIF
            IF( H .LE. YL .OR. H .GE. ZL .OR. 
     >          V .LE. YC .OR. V .GE. ZC ) THEN
              N1 = N1 + 1
C-------------- SORT = path length of a particle when stopped 
              SORT(I) = F(6,I)
              CALL KSTOP(4,I,IEX(I),*15)
 15           CONTINUE
              IF(NRES .GT. 0) THEN
                IF(IL .EQ. 2) THEN
                  IF(N1 .EQ. 1) WRITE(NRES,104)
                  CALL CNTOUR(
     >                        NOUT)
                  WRITE(NRES,106) N1,NOUT,LET(I),IEX(I)
     >            ,(FO(J,I),J=1,5),(F(J,I),J=1,5),I,IREP(I)
                  WRITE(NRES,105) F(6,I)
                ENDIF
              ENDIF
            ENDIF
          ENDIF
 
        ELSEIF(IFRM .GE. 11) THEN
C--------- Phase-space (acceptance) collimator
          GAM = (1.D0+ALP*ALP) / BET

          IF    (IFRM .EQ.11 .OR. IFRM .EQ.14) THEN
C----------- Horizontal  
            Y2 = F(2,I)*UNIT(1) - XM
            T2 = F(3,I)*UNIT(2) - XPM
            YT = Y2*T2
            Y2 = Y2*Y2
            T2 = T2*T2
            IF( GAM*Y2+2.D0*ALP*YT+BET*T2 .GT. EPSPI ) THEN
              N1 = N1 + 1
              SORT(I) = F(6,I)
              CALL KSTOP(4,I,IEX(I),*12)
 12           CONTINUE
              IF(NRES .GT. 0) THEN
                IF(IL .EQ. 2) THEN
                  IF(N1 .EQ. 1) WRITE(NRES,104)
                  CALL CNTOUR(
     >                        NOUT)
                  WRITE(NRES,106) N1,NOUT,LET(I),IEX(I)
     >            ,(FO(J,I),J=1,5),(F(J,I),J=1,5),I,IREP(I)
                  WRITE(NRES,105) F(6,I)
                ENDIF
              ENDIF
            ENDIF
 
          ELSEIF(IFRM .EQ.12 .OR. IFRM .EQ.15) THEN
C----------- Vertical  
            Y2 = F(4,I)*UNIT(3) - XM
            T2 = F(5,I)*UNIT(4) - XPM
            YT = Y2*T2
            Y2 = Y2*Y2
            T2 = T2*T2
            IF( GAM*Y2+2.D0*ALP*YT+BET*T2 .GT. EPSPI ) THEN
              N1 = N1 + 1
              SORT(I) = F(6,I)
              CALL KSTOP(4,I,IEX(I),*13)
 13           CONTINUE
              IF(NRES .GT. 0) THEN
                IF(IL .EQ. 2) THEN
                  IF(N1 .EQ. 1) WRITE(NRES,104)
                  CALL CNTOUR(
     >                        NOUT)
                  WRITE(NRES,106) N1,NOUT,LET(I),IEX(I)
     >            ,(FO(J,I),J=1,5),(F(J,I),J=1,5),I,IREP(I)
                  WRITE(NRES,105) F(6,I)
                ENDIF
              ENDIF
            ENDIF

          ELSEIF(IFRM .EQ.13 .OR. IFRM .EQ.16) THEN
C----------- Time-kineticE
CCCCC----------- Time-momentum
            Y2 = F(7,I) - XM
CCCCC            T2 = F(1,I)*UNIT(6) - XPM
            P = BORO*CL9 *F(1,I) * AMQ(2,I)
            T2 = SQRT(P*P + AMQ(1,I)*AMQ(1,I))- AMQ(1,I) - XPM
            YT = Y2*T2
            Y2 = Y2*Y2
            T2 = T2*T2
            IF( GAM*Y2+2.D0*ALP*YT+BET*T2 .GT. EPSPI ) THEN
              N1 = N1 + 1
              SORT(I) = F(6,I)
              CALL KSTOP(4,I,IEX(I),*14)
 14           CONTINUE
              IF(NRES .GT. 0) THEN
                IF(IL .EQ. 2) THEN
                  IF(N1 .EQ. 1) WRITE(NRES,104)
                  CALL CNTOUR(
     >                        NOUT)
                  WRITE(NRES,106) N1,NOUT,LET(I),IEX(I)
     >            ,(FO(J,I),J=1,5),(F(J,I),J=1,5),I,IREP(I)
                  WRITE(NRES,105) F(6,I)
                ENDIF
              ENDIF
            ENDIF

          ENDIF

        ENDIF

 1      CONTINUE
       
C        CALL REBELR(KREB3,KDUM,IMX)
C        IF(KREB3.NE.99) IMX = IMX + IMAX

        CONTINUE

        CALL CNTMXR(
     >              IMX) 
        CALL CNTOUR(
     >              NOUT) 
          IF(NRES .GT. 0 ) 
     >    WRITE(NRES,109) N1, NOUT,LIVE-N1,IMX,(100.D0*(LIVE-N1))/IMX
 109      FORMAT(/,T20,' Number of particles counted out of acceptance',
     >    /,T25,' -  at  this  collimator : ',I6,
     >    /,T25,' -  by  collimations,  from  the  beginning : ',I6,
     >    /,T25,' Overall  survival  :  ',I9,' / ',I9,'  (',G10.3,'%)')
          CALL CNTSTO(NTOT) 

      RETURN
      END
