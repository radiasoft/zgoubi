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
C  Upton, NY, 11973, USA
C  -------
      SUBROUTINE MCOBJ1(KTIR,KOUV,CINE,CENTRE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*)  KTIR(*), KOUV
      LOGICAL CINE
      DIMENSION CENTRE(*)
C     ----------------------------------------------------------
C      Sorting at random on a window (KOBJ=1) or a grid (KOBJ=2)
C      Uncorrelated coordinates
C     ---------------------------------------------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
      CHARACTER(1) LET
      INCLUDE "C.FAISCT.H"     ! COMMON/FAISCT/ LET(MXT)
      INCLUDE "C.FITEX.H"     ! COMMON/FITEX/ DN0,C0,C1,C2,C3,DL
      CHARACTER(1) KAR(41)
      INCLUDE "C.KAR.H"     ! COMMON/KAR/ KAR
      INCLUDE "C.OBJET.H"     ! COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      INCLUDE "C.REBELO.H"   ! COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.UNITS.H"     ! COMMON/UNITS/ UNIT(MXJ)

      DIMENSION AMAX(MXJ),P(MXJ),IG(MXJ),CUT(MXJ)
      DIMENSION FAISC1(MXJ,MXT),FAISC2(MXJ,MXT)

      PARAMETER (MXJ1=MXJ-1)

      DATA IKAR / 0 /

      IF(KOUV .EQ. 'Grid') THEN

C------- Read samplings of the Grid
        NL = 50
        DO 2 J= 1, MXJ1

          J2 = J+1
          IF(J2 .EQ. MXJ1+1) J2 = 1
          L2 = NL + J - 1

C------- Number of bars of the Grid :
          IG(J2) = A(NOEL,L2)

C--------- ECARTEMENT DES BARREAUX :
          P(J2) = A(NOEL,10+L2)

 2      CONTINUE

        NL=70

      ELSE

        DO 7 I=1,MXJ1
 7        IG(I) = 1

        NL=50

      ENDIF
 
C     .... LECTURE  FRONTIERES DU TIRAGE
C     .. POSITION DU BARREAU CENTRAL :
C                      Y   T   Z   P   X  D

      DO 3 J= 1, MXJ1

        J2 = J+1
        IF(J2 .EQ. MXJ1+1) J2 = 1
        L2 = NL + J - 1

C----- +/- widths (of window or bars)
        AMAX(J2) = A(NOEL,L2)

 3    CONTINUE

C------- Cut-offs
      DO 1 J = 1, MXJ1
        J2 = J+1
        IF(J2 .EQ. MXJ1+1) J2 = 1
        IF(KTIR(J) .EQ. 'Gaussian') THEN
          L2 = NL + 10 + J - 1
          CUT(J2) = A(NOEL,L2) 
        ELSE
          CUT(J2) = 1.D0
        ENDIF
 1    CONTINUE

C     .. COEFFICIENTS POUR LA LOI DE D
      DN0 = A(NOEL,NL+20)
      C0 = A(NOEL,NL+21)
      C1 = A(NOEL,NL+22)
      C2 = A(NOEL,NL+23)
      C3 = A(NOEL,NL+24)
 
C     ... LECTURE GENERATEURS ( SEULEMENT AU 1-ER PASSAGE )
C      TIRAGE AVEC MELANGE ALEATOIRE
      IF(IPASS .EQ. 1) THEN
        NI = NL+30
        IF(NI .EQ. 100) NI = 95 
        IR1 = A(NOEL,NI)
        IR2 = A(NOEL,NI+1)
        IR3 = A(NOEL,NI+2)
        IR1=(IR1/2)*2+1
        IR2=(IR2/2)*2+1
        IR3 =(IR3/2)*2+1
      ENDIF
 
      IF(NRES.GT.0) THEN
 
        WRITE(NRES,100) KOUV
100     FORMAT(15X,' Distribution in a ',A)

        IF(KOUV .EQ. 'Grid') THEN
          WRITE(NRES,101) (IG(J),J=2,MXJ1),IG(1)
 101      FORMAT(/,15X,' Number of bars of the Grid :'
     >    ,/,11X,'Y , T , Z , P , X , BR/BORO :    ',T50,6I10,/)
          WRITE(NRES,102) (P(J),J=2,MXJ1),P(1)
 102      FORMAT(15X,' Intervals between bars (MKSA units) : '
     >    ,/,11X,'Y , T , Z , P , X , BR/BORO : ',T50,1P,5E11.3,G12.4)
        ENDIF
 
        IF(CINE) WRITE(NRES,FMT=
     >  '(15X,'' Kinematical  coefficients : '',1P,G12.4
     >  ,'' /mrad'')') AMAX(1)
 
        WRITE(NRES,123) (CENTRE(J),J=2,MXJ1),CENTRE(1)
 123    FORMAT(/,15X,' Central values (MKSA units): '
     >  ,/,11X,' Yo, To, Zo, Po, Xo, BR/BORO  : ',T50,1P,6G12.4,/)
 
        WRITE(NRES,109) (AMAX(J),J=2,MXJ1) , AMAX(1)
 109    FORMAT(15X,' Width  (  +/- , MKSA units ) :'
     >  ,/,11X,' DY, DT, DZ, DP, DX, DBR/BORO : ',T50,1P,6G12.3,/)

        WRITE(NRES,110) (CUT(J),J=2,MXJ1) , CUT(1)
 110    FORMAT(15X,' Cut-offs  ( * +/-Width ) :'
     >  ,/,11X,' NY, NT, NZ, NP, NX, NBR/BORO : ',T50,6F12.2,/)

        WRITE(NRES,111) (KTIR(J),J=2,MXJ1) ,KTIR(1)
 111    FORMAT(15X,' Type  of  sorting :'
     >  ,/,11X,' Y, T, Z, P, X, D : ',T50,6A12,/)

       ENDIF
 
C--------- Change units from MKSA to cm, mrad...
       DO 33 J=1,MXJ1
          IU = J-1
          IF(IU .EQ. 0) IU=6
          CENTRE(J)=CENTRE(J)/UNIT(IU)
          P(J) = P(J)/UNIT(IU)
          AMAX(J) = AMAX(J)/UNIT(IU)
 33    CONTINUE

C------- CONSTITUTION DU FAISCEAU
C------- TIRAGE GENERATEUR IR1
       DO 31 I=1,IMAX
          DO 31 J=1,MXJ1

            CENTR=CENTRE(J)
            IF(KOUV .EQ. 'Grid') THEN
              IGR = 1 + IG(J) * RNDM() * .999999D0
              IGR = IGR/2 * (-1)**IGR
              CENTR = CENTR + IGR * P(J)
            ENDIF
 
            IF    (KTIR(J) .EQ. 'Uniform') THEN
              FAISC1(J,I)= CENTR+2.D0*(RNDM()-.5D0)*AMAX(J)

            ELSEIF(KTIR(J) .EQ. 'Gaussian') THEN
              FAISC1(J,I)=APHERF(CENTR,AMAX(J))

            ELSEIF(KTIR(J) .EQ. 'Parabolic') THEN

            ELSEIF(KTIR(J) .EQ. 'Exponentl') THEN
              FAISC1(J,I)=APHERF(CENTR,AMAX(J))

            ENDIF
 
 31     CONTINUE
 
C------- TIRAGE GENERATEUR IR2
        DO 41 I=1,IMAX
          DO 41 J=1,MXJ1
            IU = J-1
            IF(IU .EQ. 0) IU=6
            CENTR=CENTRE(J)
            IF(KOUV .EQ. 'Grid') THEN
              IGR = 1 + IG(J) * RNDM() * .999999D0
              IGR = IGR/2 * (-1)**IGR
              CENTR = CENTR + IGR * P(J)
            ENDIF
 
            IF    (KTIR(J) .EQ. 'Uniform') THEN
              FAISC2(J,I)= CENTR+2.D0*(RNDM()-.5D0)*AMAX(J)

            ELSEIF(KTIR(J) .EQ. 'Gaussian') THEN
              FAISC2(J,I)=APHERF(CENTR,AMAX(J))

            ELSEIF(KTIR(J) .EQ. 'Parabolic') THEN

            ELSEIF(KTIR(J) .EQ. 'Exponentl') THEN
              FAISC2(J,I)=APHERF(CENTR,AMAX(J))

            ENDIF
 
 41     CONTINUE
 
      CALL MIX(IMAX,FO,FAISC1,FAISC2)
 
      XL = 0.D0
      DO 5 I=1,IMAX
        IF(CINE) FO(1,I)=CENTRE(1) + FO(3,I)*AMAX(1)
        DO 4 J=1,MXJ1
          F(J,I)=FO(J,I)
 4      CONTINUE
        IKAR = IKAR + 1
        IF(IKAR.GT.41) IKAR=1
        LET(I)= KAR(IKAR)
        IREP(I)=I
        IEX(I)=1
        XX = ABS( F(6,I) )
        IF(XL .LT. XX ) XL = XX
 5    CONTINUE
 
C     ... Long target : projct coordinates 
C         onto plane orthogonal to longitudinal axis 
C         at centre of target
      IF ( XL .NE. 0.D0) THEN
        DO 6 I=1,IMAX
           XX = F(6,I)
           F(2,I)=F(2,I) - XX*TAN(F(3,I)*1.D-3)
           F(4,I)=F(4,I) - XX*TAN(F(5,I)*1.D-3)/COS(F(3,I)*1.D-3)
           F(6,I)= XX/(COS(F(3,I)*1.D-3)*COS(F(5,I)*1.D-3))
C           F(6,I)= - AL/(COS(F(3,I)*1.D-3)*COS(F(5,I)*1.D-3))
 6      CONTINUE
      ENDIF
 
      RETURN
      END
