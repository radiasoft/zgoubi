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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE CIBLE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CINE/ M1,M2,M3,M4,M12,M1212,M22,P0 ,G,C,C2,BG,EL3M,PC3
     1,THETA,BETA,Q,PS(5),TS(5),NPS,NTS,II
      COMMON/CONST2/ ZERO, UN
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      CHARACTER LET
      COMMON/FAISCT/ LET(MXT)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      COMMON/DON/ A(09876,99),IQ(09876),IP(09876),NB,NOEL
      COMMON/OBJET/ FO(MXJ,MXT),KOBJ,IDMAX,IMAXT
      LOGICAL ZSYM
      COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT
 
      DOUBLE PRECISION M1,M2,M3,M4,M12,M1212,M22
 
C         M1  MASSE CIBLE EN MEV
C           M2 MASSE PROJECTILE EN MEV
C        M3 MASSE DU DIFFUSE EN MEV
C           M4    RESIDUEL  EN MEV
C        Q   Q REACTION
C            T2     ENERGIE DU PROJECTILE
C           THETA ANGLE DE DIFFUSION
C          BETA ANGLE CIBLE
C        LES ANGLES SONT DONNES EN DEGRES
C         LE PGM CONVERTIT LES DONNES EN GEV
 
      M1 = A(NOEL,1)
      M2 = A(NOEL,2)
      M3 = A(NOEL,3)
      Q  = A(NOEL,4)
      T2 = A(NOEL,5)
      THETA = A(NOEL,6)
      BETA = A(NOEL,7)
 
       M1=M1/1.D3
      M2=M2/1.D3
      M3=M3/1.D3
       T2=T2/1.D3
       P0=SQRT(T2*T2+2.D0*T2*M2)
      IF(NRES .GT. 0) WRITE(NRES,100)
 100   FORMAT('1',10('*'),' CIBLE ',/)
      IF(NRES .GT. 0) WRITE(NRES,1000) M1,M2,M3,Q,T2
 1000  FORMAT(5X,' DIFFUSION ',/,5X,' M1= ',G15.8,' GEV',3X
     1,' M2= ',G15.8,' GEV',3X,' M3= ',G15.8,' GEV',/,5X,' Q = ',
     2G15.8,3X,' T2= ',G15.8)
      IF(NRES .GT. 0) WRITE(NRES,1001) THETA,BETA
 1001  FORMAT(/,5X,' ANGLE DE DIFFUSION ',G15.8,' DEG ',3X,
     1' ANGLE CIBLE  ',G15.8,' DEG')
      THETA=THETA*0.0174533D0
      BETA=BETA*0.0174533D0
      M12=M1+M2
      M1212=M12*M12
      M22=M2*M2
 
      NTSMAX = A(NOEL,10)
      NPSMAX = A(NOEL,11)
      PTS = A(NOEL,20)
      PPS = A(NOEL,21)
      DTS = A(NOEL,22)
 
      DTS=DTS*1.D-3
      IF(NRES .GT. 0) WRITE(NRES,101) NTSMAX,PTS,NPSMAX,PPS
      PTS=PTS*1.D-3
      PPS=PPS*1.D-3
      DO 2 NTS=1,NTSMAX
        TS(NTS)=DELTA(NTS,PTS)
        TS(NTS)=TS(NTS)+DTS
 2    CONTINUE
      DO 3 NPS=1,NPSMAX
        PS(NPS)=DELTA(NPS,PPS)
 3    CONTINUE
      NSMAX=NTSMAX*NPSMAX
      DO 4 I=1,IMAX
        II=IMAX-I+1
        IF=1+(II-1)*NSMAX
        IEX(IF)=IEX(II)
        F(6,IF)= F(6,II)
        LET(IF)=LET(II)
        IIREP=IREP(II)
        IREP(IF)=1+(IIREP-1)*NSMAX
        DO 18 J=1,5
          FO(J,IF)=FO(J,II)
          F(J,IF)=F(J,II)
   18   CONTINUE
        DO 5 NS=2,NSMAX
          IF=IF+1
          F(6,IF)= F(6,II)
          DO 19 J=1,5
            FO(J,IF)=FO(J,II)
   19     CONTINUE
          LET(IF)=LET(II)
5       CONTINUE
4     CONTINUE
      F(6,IF)= F(6,II)
      IMAX=IMAX*NSMAX
      IMAXT=IMAXT*NSMAX
      DO 1 ID=1,IDMAX
        IMAX1=1+(ID-1)*IMAXT
        IMAX2=IMAX1+IMAXT-NSMAX
        DP=F(1,IMAX1)
        CALL INICIN
        DO 6 I=IMAX1,IMAX2,NSMAX
          II=I-1
          IF(IREP(I) .NE. I)  GOTO 12
          CALL INITRA(I)
          CALL CHAREF(.FALSE.,ZERO,ZERO,BETA)
          PO=P
          Y0=Y
          T0=T
          Z0=Z
          IF(F(5,I) .EQ. 0.D0 .AND. F(4,I) .EQ. 0.D0)  GOTO 9
          DO 7 NTS=1,NTSMAX
            DO 8 NPS=1,NPSMAX
              Z=Z0
              Y=Y0
              P=PO
              T=T0
              II=II+1
              IREP(II)=II
              CALL CINEMA
8           CONTINUE
7         CONTINUE
          GOTO 15
9         DO 10 NTS=1,NTSMAX
            DO 11 NPS=1,NPSMAX,2
              Z=Z0
              Y=Y0
              P=PO
              T=T0
              II=II+1
              IREP(II)=II
              CALL CINEMA
              IF(NPSMAX .EQ. 1) GOTO 10
              IF(NPS .EQ. 1)  GOTO 11
              II=II+1
              IREP(II)=II-1
              CALL DEJACA(II)
11          CONTINUE
10        CONTINUE
          GOTO 15
12        CONTINUE
          IR=IREP(I)-1
          DO 13 NTS=1,NTSMAX
            DO 14 NPS=1,NPSMAX
              II=II+1
              IR=IR+1
              IREP(II)=IR
              IF(PS(NPS).GT.0.D0)  IREP(II)=IR+1
              IF(PS(NPS).LT.0.D0)  IREP(II)=IR-1
              CALL DEJACA(II)
 14         CONTINUE
 13       CONTINUE
 
 15       CONTINUE
 
6       CONTINUE
1     CONTINUE
 
C        AJUSTEMENT LIGNE EN AVAL DE LA CIBLE
C       ***************************************************************
 
       IF(NRES .GT. 0) WRITE(NRES,196)
 196   FORMAT(//,5X,' AJUSTEMENT DE LA LIGNE AVAL CIBLE :',/)
 
       BORO = A(NOEL,30)
 
       IF(NRES .GT. 0) WRITE(NRES,197) BORO
 197   FORMAT(10X,'MAGNETIC RIGIDITY AFTER TARGET',3X,G15.8,2X,
     >  'kG*cm')
       DO 199 I=2,IMAX
       F(1,I)=F(1,I)/F(1,1)
 199   CONTINUE
       F(1,1)=1D0
      RETURN
  101  FORMAT(/,5X,I4,' ANGLES EN THETA , PAS = ',G8.1,' MRAD',/
     > ,5X,I5,' ANGLES EN PHI , PAS = ',G8.1,' MRAD',/)
      END
