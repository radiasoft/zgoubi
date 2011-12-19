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
C  Brookhaven National Laboratory               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE OBJ6
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ***************************************
C     CONSTITUTION DE L'OBJET INITIAL KOBJ=60
C     ***************************************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IIP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,40)
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
C----- KAR: LETTRES AFFECTEES AUX TRAJECTOIRES ( 'S'  EST RESERVEE
C      POUR ETIQUETER LES PARTICULES SECONDAIRES -OPTION 'MCDESINT')
      CHARACTER  KAR(41)
      COMMON/KAR/ KAR
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
 
      PARAMETER(MXJ1=MXJ-1)

      DIMENSION IDE(5),JDE(5),P(MXJ)
      EQUIVALENCE (IDE(2),IYMAX),(IDE(3),ITMAX),(IDE(4),IZMAX),
     > (IDE(5),IPMAX),(IDE(1),IMAXD)
      EQUIVALENCE (JDE(2),IY   ),(JDE(3),IT   ),(JDE(4),IZ   ),
     > (JDE(5),IP   ),(JDE(1),ID)
 
      DIMENSION REF(MXJ)

      P(2) = A(NOEL,20)
      P(3) = A(NOEL,21)
      P(4) = A(NOEL,22)
      P(5) = A(NOEL,23)
      P(6) = A(NOEL,24)
      P(1) = A(NOEL,25)
      IDMAX= 1

      K = 30
      DO 62 J = 2,MXJ1
        REF(J) = A(NOEL,K)
        K = K + 1
 62   CONTINUE
      REF(1) = A(NOEL,35)

C     ** TRAJECTOIRES SANS COUPLAGE
      II=1
      DO 69 J=2,5
        DX = P(J)
        DO 69 I=1,4
          II = II+1
          I2 = I/2
          IF(2*I2 .EQ. I) THEN
            FO(J,II) =  - I2 * DX
          ELSE
            FO(J,II) =  + (1+I2) * DX
          ENDIF
 69   CONTINUE

      DO 67 I=1,4
        II = II+1
        I2 = I/2
        IF(2*I2 .EQ. I) THEN
          FO(1,II) =  - I2 * P(1)
        ELSE
          FO(1,II) =  + (1+I/2) * P(1)
        ENDIF
 67   CONTINUE
 
C     ** TRAJECTOIRES AVEC COUPLAGE ./YT, /YZ, /YP, YD
      J=2
      DO 61 JJ=3,5
        II=II+1
        FO(J,II) =  + P(J)
        FO(JJ,II) =  + P(JJ)
        II=II+1
        FO(J,II) =  + P(J)
        FO(JJ,II) =  - P(JJ)
        II=II+1
        FO(J,II) =  - P(J)
        FO(JJ,II) =  - P(JJ)
        II=II+1
        FO(J,II) =  - P(J)
        FO(JJ,II) =  + P(JJ)
 61   CONTINUE
      II=II+1
      FO(J,II) =  + P(J)
      FO(1,II) =  + P(1)
      II=II+1
      FO(J,II) =  + P(J)
      FO(1,II) =  - P(1)
      II=II+1
      FO(J,II) =  - P(J)
      FO(1,II) =  - P(1)
      II=II+1
      FO(J,II) =  - P(J)
      FO(1,II) =  + P(1)
 
C     ** TRAJECTOIRES AVEC COUPLAGE ./TZ, ./TP ET ./TD
      J=3
      DO 63 JJ=4,5
        II=II+1
        FO(J,II) =  + P(J)
        FO(JJ,II) =  + P(JJ)
        II=II+1
        FO(J,II) =  + P(J)
        FO(JJ,II) =  - P(JJ)
        II=II+1
        FO(J,II) =  - P(J)
        FO(JJ,II) =  - P(JJ)
        II=II+1
        FO(J,II) =  - P(J)
        FO(JJ,II) =  + P(JJ)
 63   CONTINUE
      II=II+1
      FO(J,II) =  + P(J)
      FO(1,II) =  + P(1)
      II=II+1
      FO(J,II) =  + P(J)
      FO(1,II) =  - P(1)
      II=II+1
      FO(J,II) =  - P(J)
      FO(1,II) =  - P(1)
      II=II+1
      FO(J,II) =  - P(J)
      FO(1,II) =  + P(1)
 
C     ** TRAJECTOIRES AVEC COUPLAGE ./ZP ET ./ZD
      J=4
      JJ=5
      II=II+1
      FO(J,II) =  + P(J)
      FO(JJ,II) =  + P(JJ)
      II=II+1
      FO(J,II) =  + P(J)
      FO(JJ,II) =  - P(JJ)
      II=II+1
      FO(J,II) =  - P(J)
      FO(JJ,II) =  - P(JJ)
      II=II+1
      FO(J,II) =  - P(J)
      FO(JJ,II) =  + P(JJ)
      II=II+1
      FO(J,II) =  + P(J)
      FO(J,II) =  + P(J)
      FO(1,II) =  + P(1)
      II=II+1
      FO(J,II) =  + P(J)
      FO(1,II) =  - P(1)
      II=II+1
      FO(J,II) =  - P(J)
      FO(1,II) =  - P(1)
      II=II+1
      FO(J,II) =  - P(J)
      FO(1,II) =  + P(1)
 
C     ** TRAJECTOIRES AVEC COUPLAGE  ./PD
      II=II+1
      FO(5,II) =  + P(5)
      FO(1,II) =  + P(1)
      II=II+1
      FO(5,II) =  + P(5)
      FO(1,II) =  - P(1)
      II=II+1
      FO(5,II) =  - P(5)
      FO(1,II) =  - P(1)
      II=II+1
      FO(5,II) =  - P(5)
      FO(1,II) =  + P(1)
 
      IMAX = II
      IMAXT=IMAX/IDMAX
      IREP(1)=1
      IKAR=0
      DO 66 I=1,IMAX
         IEX(I) =1
         IREP(I)=I
         IKAR = IKAR+1
         IF(IKAR.GT.41)  IKAR=1
         LET(I)=KAR(IKAR)
         DO 66 J=1,6
            FO(J,I) = FO(J,I) + REF(J)
            F(J,I)=FO(J,I)
 66   CONTINUE
      IF(NRES.GT.0) THEN
        WRITE(NRES,100) KOBJ,IMAX
        WRITE(NRES,102) (P(J),J=1,MXJ1)
        WRITE(NRES,FMT='(/,19X,'' Reference trajectory : '',
     >    4X,F7.4,6(4X,F6.2),/)') (REF(I), I=1, MXJ1)
      ENDIF
  100 FORMAT(/,41X,'CALCUL  DES  TRAJECTOIRES',//,30X,'OBJET  (',I1,
     1')  FORME  DE ',I6,' POINTS ',//)
  102 FORMAT(/,19X,' ECHANTILLONNAGE ',4X,F7.4,5(4X,F6.2),/)
      RETURN
      END
