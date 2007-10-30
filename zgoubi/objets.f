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
      SUBROUTINE OBJETS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ****************************************
C     CONSTITUTION DE L'OBJET INITIAL STOCKE DANS F0(6,MXT):
C       1-ERE INDICE DE F0: D, Y0, T0, Z0, P0 POUR CHAQUE TRAJECTOIRE,
C       2-EME INDICE DE F0: NUMERO DE TRAJECTOIRE .
C     EN COURS DE VOL, F(MXJ,MXT) STOCKE LES COORDONNEES CALCULEES
C     ****************************************
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
C      INCLUDE "MAXTRA.H"
C      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
C     > ,YCH,ZCH
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IIP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,20)
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT)
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
C----- KAR: tagging letter ( 'S'  is reserved for tagging secondary particles 
C            as resulting from decay (keyword 'MCDESINT')
      CHARACTER KAR(41)
      COMMON/KAR/ KAR
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      COMMON/RIGID/ BORO,DPREF,DP,BR
      COMMON/SYNCH/ RET(MXT), DPR(MXT),PS
 
      PARAMETER(MXJ1=MXJ-1)
      PARAMETER (NTM=41)
      DIMENSION DE(5,NTM),IDE(5),JDE(5),P(MXJ)
      EQUIVALENCE (IDE(2),IYMAX),(IDE(3),ITMAX),(IDE(4),IZMAX),
     > (IDE(5),IPMAX),(IDE(1),IMAXD)
      EQUIVALENCE (JDE(2),IY   ),(JDE(3),IT   ),(JDE(4),IZ   ),
     > (JDE(5),IP   ),(JDE(1),ID)
 
      DIMENSION REF(MXJ)

C----- MAGNETIC  RIGIDITY (KG*CM), MASS (MeV/c/2)
      BORO = A(NOEL,1)
      IF(NRES.GT.0) WRITE(NRES,103) BORO
 103  FORMAT(25X,' MAGNETIC  RIGIDITY =',F15.3,' kG*cm')

C----- Get data from possible REBELOTE keyword, then some action
      CALL REBELR(KREB3,KREB31)
      IF(KREB3 .EQ. 99) THEN
C------- Set to 99 by REBELOTE
C------- For simulation of multiturn injection
        IF(NRES.GT.0) WRITE(NRES,133) IMAX
 133    FORMAT(//,15X,'Final  coordinates  of  previous  run',1X
     >  ,' taken  as  initial  coordinates ; ',I9,' particles')
        IF(KREB31 .NE. 0) THEN
C--------- add new beamlet next to the previous one(s), e.g. for multiturn injection
          IF(IPASS .LE. 1+KREB31) THEN
            IF(NRES.GT.0) WRITE(NRES,FMT='(
     >           15X,'' Injection run ; new beamlet launched'',/)')
          ELSE
            GOTO 99
          ENDIF
        ELSE
          GOTO 99
        ENDIF
      ENDIF

C      CALL RAZ(FO,MXJ*MXT)
C----- Was necessary for FIT procedure when time is constrained : 
      CALL RAZ(F,MXJ*MXT)

      KOBJ = A(NOEL,10)
      KOBJ2 = NINT(10.D0*A(NOEL,10)) - 10*KOBJ

      GOTO( 1, 2,16,97,50,60,1 ,80),KOBJ

      IF(NRES.GT.0) WRITE(NRES,FMT=
     >'(10X,''SBR OBJETS:    NO  SUCH  OBJECT  KOBJ='',I2)') KOBJ
      CALL ENDJOB(' NO  SUCH  OBJECT  ',-99)
 
C---------- OBJET with 11 traj. for 1st order matrix calculation
 50   CONTINUE
      CALL RAZ(FO,MXJ*MXT)
C      CALL OBJ52(KOBJ2)
      CALL OBJ5
      CALL MATRI2(KOBJ2)
      GOTO 99

C---------- OBJET with 61 traj. for 1st, 2nd, and higher coeff. computation
 60   CONTINUE
      CALL RAZ(FO,MXJ*MXT)
      CALL OBJ6
      GOTO 99

C---------- Read OBJET from file
 16   CONTINUE
        CALL OBJ3(KOBJ2,BORO)
      GOTO 99

C---------- Initial conditions particle by particle
 2    CONTINUE
      IMAX = A(NOEL,20)
      IF(IMAX .GT. MXT) GOTO 98
      IDMAX= A(NOEL,21)

C----- For allowing the use of the first 7 traj with FIT
      II = 20
      DO 11 I=1,7
        II = II + 10
        FO(1,I) = A(NOEL,II+5)
        DO 112 J=2,MXJ1
          FO(J,I) = A(NOEL,II+J-2)
 112      CONTINUE
C        P0 = BORO*CL9*FO(1,I)
C        BTA = P0 / SQRT( P0*P0 + 0.511**2 )  
C        FO(7,I) = FO(6,I)/(BTA*CL)
 11   CONTINUE

      LUN = NRES
      IF(NRES.LE.0) LUN=6
 
      DO 3 I=1,IMAX
        IF(FO(1,I) .EQ. 0.D0) THEN
          IEX(I) = -6
            WRITE(LUN,*) ' Momentum value 0 found, ', 
     >          ' particle of concern  is  # ',I,' ; ', 
     >          ' its KEX will be forced to -6 ; will not be tracked'
            WRITE(LUN,*) 'Y T Z P S D Time : ',(F(J,I),J=1,7),
     >                                            ' KEX=',IEX(I)
        ENDIF
C Time of flight is initialized in subroutine PARTIC
C        P0 = BORO*CL9*FO(1,I)
C        BTA = P0 / SQRT( P0*P0 + 0.511**2 )  
C        FO(7,I) = FO(6,I)/(BTA*CL)
        DO 4 J=1,MXJ
          F(J,I)=FO(J,I)
 4      CONTINUE
        IREP(I)=I
    3 CONTINUE

      IMAXT=IMAX/IDMAX
      IF(NRES.GT.0) WRITE(NRES,106) KOBJ,IMAX
      GOTO 99
 
C---------- Initial conditions on an ellipsoid
 80   CONTINUE
       CALL OBJ8(KREB3,KREB31)
      GOTO 99
 

C---------- OBJET AUTOMATIQUE SYMETRIQUE
 1    CONTINUE
      IYMAX = A(NOEL,20)
      ITMAX = A(NOEL,21)
      IZMAX = A(NOEL,22)
      IPMAX = A(NOEL,23)
      IXMAX = A(NOEL,24)
      IDMAX = A(NOEL,25)
      IF(    IYMAX .GT. NTM .OR. ITMAX .GT. NTM .OR. IZMAX .GT. NTM
     > .OR. IPMAX .GT. NTM .OR. IXMAX .GT. NTM .OR. IDMAX .GT. NTM) THEN
        CALL ENDJOB('Too many trajectories, max is ',MXT)
      ENDIF
      IMAXD=IDMAX
      P(2) = A(NOEL,30)
      P(3) = A(NOEL,31)
      P(4) = A(NOEL,32)
      P(5) = A(NOEL,33)
      P(6) = A(NOEL,34)
      P(1) = A(NOEL,35)
      D = A(NOEL,45)

      IF(KOBJ2.EQ.0) THEN
      DO 6 J=1,5
         IDEMAX=IDE(J)
         DO 7 K=1,IDEMAX
C           DE(J,K)=DELTA(K,P(J))
            K2=K/2
            IF(2*K2 .EQ. K) THEN
               DE(J,K) = K2*P(J)
            ELSE
               DE(J,K) =-K2*P(J)
            ENDIF
    7    CONTINUE
    6 CONTINUE
      DO 14 K=1,IDMAX
         DE(1,K)=DE(1,K)+D
   14 CONTINUE
      I=0
      DO  8 ID=1,IDMAX
        IKAR=0
        DO  8 IY=1,IYMAX
          DO  8 IT=1,ITMAX
            IKAR=IKAR+1
            IF(IKAR.GT.41)  IKAR=1
            DO  8 IZ=1,IZMAX
              DO  8 IP=1,IPMAX
                I=I+1
                IREP(I)=I
                IF(IZ .EQ. 1 .AND. DE(5,IP).LT.0.D0)  IREP(I)=I-1
                IF(DE(4,IZ).LT.0D0 .AND. DE(5,IP) .EQ. 0.D0)
     >           IREP(I)=I-IPMAX
                IF(DE(4,IZ).LT.0D0 .AND. DE(5,IP) .GT. 0.D0)
     >           IREP(I)=I-IPMAX+1
                IF(DE(4,IZ).LT.0D0 .AND. DE(5,IP) .LT. 0.D0)
     >           IREP(I)=I-IPMAX-1
                DO 13 J=1,5
                  KDE=JDE(J)
                  FO(J,I)=DE(J,KDE)
C                  F(J,I)=FO(J,I)
                  LET(I)=KAR(IKAR)
13              CONTINUE
        F(6,I)= 0D0
8     CONTINUE
      IMAX=I
 
      IF(IMAX .GT. MXT) GOTO 98
 
      K = 40
      DO 162 J = 2,MXJ1
        REF(J) = A(NOEL,K)
        K = K + 1
 162  CONTINUE

      REF(1) = 0.D0
      DO 161 I=1,IMAX
        DO 161 J=1, 6
          FO(J,I) = FO(J,I) + REF(J)
          F(J,I)=FO(J,I)
 161  CONTINUE

      ELSEIF(KOBJ2.EQ.1) THEN
C--------------- OBJET with Z>0 ET P>0
      DO 21 J=1,3
         IDEMAX=IDE(J)
         DO 22 K=1,IDEMAX
            K2=K/2
            IF(2*K2 .EQ. K) THEN
               DE(J,K) = K2*P(J)
            ELSE
               DE(J,K) =-K2*P(J)
            ENDIF
22       CONTINUE
21    CONTINUE
      DO 23 J=4,5
         IDEMAX=IDE(J)
         DO 23 K=1,IDEMAX
            DE(J,K) =(K-1)*P(J)
23    CONTINUE
      DO 29 K=1,IDMAX
29       DE(1,K)=DE(1,K)+D
      I=0
      DO 25 ID=1,IDMAX
         IKAR=0
         DO 25 IY=1,IYMAX
            DO 25 IT=1,ITMAX
               IKAR=IKAR+1
               IF(IKAR.GT.41)  IKAR=1
               DO 25 IZ=1,IZMAX
                  DO 25 IP=1,IPMAX
                     I=I+1
                     IREP(I)=I
                     DO 26 J=1,5
                        KDE=JDE(J)
                        FO(J,I)=DE(J,KDE)
                        F(J,I)=FO(J,I)
                        LET(I)=KAR(IKAR)
26                   CONTINUE
                  F(6,I)= 0D0
25    CONTINUE
      IMAX=I

      IF(IMAX .GT. MXT) GOTO 98

      K = 40
      DO 24 J = 2,MXJ1
        REF(J) = A(NOEL,K)
        K = K + 1
 24   CONTINUE

      REF(1) = 0.D0
      DO 28 I=1,IMAX
        DO 28 J=1, 6
          FO(J,I) = FO(J,I) + REF(J)
          F(J,I)=FO(J,I)
 28   CONTINUE

      IMAXT=IMAX/IDMAX
      IF(NRES.GT.0) THEN
        WRITE(NRES,106) KOBJ,IMAX
        WRITE(NRES,101) (IDE(J),J=1,MXJ1)
        WRITE(NRES,102) (P(J),J=1,MXJ1)
      ENDIF
      DO 27 I=1,IMAX
         IEX(I) = 1
   27 CONTINUE
      ENDIF


      IF(KOBJ .EQ. 7) THEN
C------- EFFET CINEMATIQUE PRIS EN COMPTE
        DO 19 I=1,IMAX
          FO(1,I)=D + FO(3,I)*P(1)
 19       F(1,I)=FO(1,I)
      ENDIF
 
      IMAXT=IMAX/IDMAX
      IF(NRES.GT.0) THEN
        WRITE(NRES,106) KOBJ,IMAX
        WRITE(NRES,101) (IDE(J),J=1,MXJ1)
        WRITE(NRES,102) (P(J),J=1,MXJ1)
      ENDIF
      DO 15 I=1,IMAX
        IEX(I) = 1
   15 CONTINUE
      GOTO 99
 
 97   CONTINUE
      CALL ENDJOB(' NO SUCH OBJET KOBJ = ',4)

 98   CONTINUE
      CALL OBJERR(ABS(NRES),2,MXT,'    Too many particles.')
      CALL ENDJOB(' Too many particles ',-99)

 99   CONTINUE
      IF(IPASS.EQ.1) CALL CNTMXW(IMAX)
      LUN = NRES
      IF(NRES.LE.0) LUN=6
      DO 991 I=1,IMAX
        IF(IEX(I) .LT. -1) THEN
          CALL KSTOP(ABS(IEX(I)),IT,IEX(I),*992)
 992      CONTINUE
        ELSEIF(F(1,I) .EQ. 0.D0) THEN
          IEX(I) = -6
          CALL OBJERR(ABS(NRES),1,MXT,'   momentum value 0 found.')
          WRITE(LUN,*) ' particle of concern  is  # ',I,' ; ', 
     >          ' its KEX will be forced to -6 ; will not be tracked'
          WRITE(LUN,*) 'Y T Z P S D Time : ',(F(J,I),J=1,7),
     >                                            ' KEX=',IEX(I)
        ENDIF
 991  CONTINUE   
      RETURN

 
  106 FORMAT(/,41X,'CALCUL  DES  TRAJECTOIRES',//,30X,'OBJET  (',I1,
     1')  FORME  DE ',I6,' POINTS ',//)
  101 FORMAT(/,42X,'D',7X,'Y(cm)',5X,'T(mrd)',4X,'Z(cm)',5X,'P(mrd)',4X,
     >'X(cm)',//,30X,'NOMBRE',5(5X,I3,2X),/)
  102 FORMAT(/,19X,' ECHANTILLONNAGE ',4X,F6.4,5(4X,F6.2),/)
      END
