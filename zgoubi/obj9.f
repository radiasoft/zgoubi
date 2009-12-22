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
      SUBROUTINE OBJ9
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------------------
C      CONSTITUTION DE L'OBJET SUR UN ELLIPSOID
C     -----------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      COMMON/FITEX/ DN0,C0,C1,C2,C3,DL
      CHARACTER  KAR(41)
      COMMON/KAR/ KAR
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IMAXD,IMAXT
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/UNITS/ UNIT(MXJ)

      DIMENSION ALP(MXJ),BET(MXJ),EPS(MXJ)
      DIMENSION CENTRE(MXJ)
      PARAMETER(MXJ1=MXJ-1)
      CHARACTER*80 TXT
      DIMENSION ISAM(3)
      PARAMETER (MRHO=4)
      DIMENSION RHO(MRHO)

      DATA RHO / 0.2671D0, 0.9400D0, 1.9617D0, 4.1589D0 / 

C----- IX, IY, IZ
      IMAX = 0
      DO I= 1, 3
        ISAM(I) = A(NOEL,20+I-1)
        IF(ISAM(I).NE.0) IMAX = IMAX + ISAM(I) * MRHO
      ENDDO

      IF(IMAX .GT. MXT) GOTO 98
      IMI  = 1
      IMA = IMAX

C----- CENTRE ELLIPSES
C       2-Y, 3-T, 4-Z, 5-Z, 6-X, 1-D
      CENTRE(2) = A(NOEL,30)
      CENTRE(3) = A(NOEL,31)
      CENTRE(4) = A(NOEL,32)
      CENTRE(5) = A(NOEL,33)
      CENTRE(6) = A(NOEL,34)
      CENTRE(1) = A(NOEL,35)
C----- PARAMETRE ELLIPSES
      ALP(2)=A(NOEL,40)
      ALP(4)=A(NOEL,50)
      ALP(6)=A(NOEL,60)
      BET(2)=A(NOEL,41)
      BET(4)=A(NOEL,51)
      BET(6)=A(NOEL,61)
      EPS(2)=A(NOEL,42)
      EPS(4)=A(NOEL,52)
      EPS(6)=A(NOEL,62)
 
      IF(NRES.GT.0) THEN
        WRITE(NRES,100) (rho(i),i=1,4)
100     FORMAT(/,15X,' Distribution  on  four  ellipses  at  ',
     >  '  rho/epsilon = ',1P,4(1X,e12.4),/) 
 
        WRITE(NRES,123) (CENTRE(J),J=2,MXJ1),CENTRE(1)
 123    FORMAT(15X,'  Ellipse centres (m-rad): '
     >  ,/,11X,' HORIZONTAL    ( Yo, To ):',T50,1P,2G12.4
     >  ,/,11X,' VERTICAL      ( Zo, Po ):',T50,   2G12.4
     >  ,/,11X,' LONGITUDINAL  ( Xo, Do ):',T50,   2G12.4,/)
 
        WRITE(NRES,109) (ALP(J),BET(J),EPS(J),J=2,MXJ1,2)
 109    FORMAT(15X,
     >  ' Alpha(rad), Beta(m/rad), E/pi(m.rad) :'
     >  ,/,11X,' HORIZONTAL    :',T35,1P,3G12.4
     >  ,/,11X,' VERTICAL      :',T35,   3G12.4
     >  ,/,11X,' LONGITUDINAL  :',T35,   3G12.4,/)

        WRITE(NRES,FMT='(15X,'' Sampling of phases in X/Y/Z : '', 
     >             3(I5,''/''))') (ISAM(I),I=1,3)
      ENDIF

C----- CONSTITUTION DU FAISCEAU
      IT = 0
      DO 1 J=2,MXJ1,2
C For each coordinate pair ((x,x'), ...)
        J1=J+1
        J2=J/2
        IF(J1.EQ.MXJ) J1=1 
        IF(EPS(J).EQ.0.D0) THEN
          DO 11 I=IMI,IMA
            FO(J ,I)=0.D0
            FO(J1,I)=0.D0
 11       CONTINUE
        ELSE
          IF(ISAM(J2) .GT. 0) THEN 
C If sample index is not zero
            DA = 2.D0*PI/DBLE(ISAM(J2))
            DO II = 1, 4
              REB=RHO(II)*SQRT(EPS(J)*BET(J))
              ANG = 0.D0
              DO I=1,ISAM(J2)
                IT = IT+1
                X = REB*COS(ANG)
                FO(J ,IT) = X/UNIT(J-1)
                FO(J1,IT) = (REB*SIN(ANG)-ALP(J)*X)/BET(J)/UNIT(J)
              seb = sqrt(bet(j))
       xn = FO(J,IT)/seb*UNIT(J-1)
       xpn = (alp(j)*FO(J,IT)*UNIT(J-1) + bet(j)*FO(J1,IT)*UNIT(J))/seb
           write(88,*) xn, xpn
     >         ,xn*xn+xpn*xpn
C     >  ,(1+alp(j)*alp(j))/bet(j)*x*x+2.d0*alp(j)*x*xp+bet(j)*xp*xp
     > ,ang, reb,i,it
                ANG = ANG + DA 
              ENDDO
            ENDDO
          ELSE
            DO I=1,ISAM(J2)
              FO(J ,I) = 0.D0
              FO(J1,I) = 0.D0
            ENDDO
          ENDIF
        ENDIF
 1    CONTINUE      
 
      IKAR = 1
      DO 5 I=IMI,IMA
        FO(1,I)=FO(1,I) + CENTRE(1)/UNIT(6)
        F(1,I)=FO(1,I)
        DO 4 J=2,MXJ1
          FO(J,I)=FO(J,I) + CENTRE(J)/UNIT(J-1)
          F(J,I)=FO(J,I)
 4      CONTINUE
        IREP(I)=I
        IEX(I)=1
        LET(I)= KAR(IKAR)
        IKAR=IKAR+1
        IF(IKAR.GT.41) IKAR=1
 5    CONTINUE
 
      RETURN

 98   TXT = '    Too  many  particles'
      CALL OBJERR(ABS(NRES),2,MXT,TXT)
      CALL ENDJOB('OBJERR,  Too  many  particles ',-99)

      RETURN
      END
