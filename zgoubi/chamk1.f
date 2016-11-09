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
      SUBROUTINE CHAMK1(A1,R1,Z1,IMAP,*)
      USE DYNHC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ---------------------------------------------------------
C     CALCUL DU Champ ET DE SES DERIVEES, A PARTIR D'UNE CARTE,
C     AU POINT COURANT.
C      ID=1 (ZGOUBI VERSION CARTES 2-DIM)  OU 3 (CARTES 3-DIM)
C      SI ID=1 ALORS IZ=1 ;
C      SI ID=3 ALORS IZ >= NOMBRE DE CARTES EN Z= IMPAIR > 2
C     --------------------------------------------------------
      INCLUDE 'PARIZ.H'
      INCLUDE "XYZHC.H"
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH
 
      INCLUDE "C.CHAMP.H"     ! COMMON/CHAMP/ BZ0(5,5), EZ0(5,5)
      INCLUDE "C.CHAVE_2.H"     ! COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      INCLUDE "C.DDBXYZ.H"     ! COMMON/DDBXYZ/ DB(3,3),DDB(3,3,3)
      INCLUDE "C.D3B_2.H"     ! COMMON/D3BXYZ/ D3BX(3,3,3), D3BY(3,3,3), D3BZ(3,3,3)
      INCLUDE "C.D4B.H"     ! COMMON/D4BXYZ/ D4BX(3,3,3,3) ,D4BY(3,3,3,3) ,D4BZ(3,3,3,3)
      INCLUDE "C.DDEXYZ_2.H"     ! COMMON/DDEXYZ/ DE(9),DDE(27)
      INCLUDE "C.D3E_2.H"     ! COMMON/D3EXYZ/ D3EX(3,3,3), D3EY(3,3,3), D3EZ(3,3,3)
      INCLUDE "C.D4EXYZ.H"     ! COMMON/D4EXYZ/ D4EX(3,3,3,3) ,D4EY(3,3,3,3) ,D4EZ(3,3,3,3)
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      INCLUDE "C.MARK.H"     ! COMMON/MARK/ KART,KALC,KERK,KUASEX
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.ORDRES.H"     ! COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
 
      DIMENSION A000(3),A100(3),A010(3),A001(3),A200(3),A020(3),
     >A002(3),A110(3),A101(3),A011(3)
      PARAMETER (F1S6=1.D0/6.D0, F1S12=1.D0/12.D0, F1S24=1.D0/24.D0,
     >F1S35=1.D0/35.D0, F1S49=1.D0/49.D0, F1S70=1.D0/70.D0)
 
      DIMENSION BC(2),DBC(2,2),DDBC(2,2,2)
 
      PARAMETER (MDX=6)
      DIMENSION BX(MDX)

      DIMENSION BMESH(5,5),BMESH3(3,3,3),BMESH1(9)
 
      EQUIVALENCE (BMESH1(1),F1),(BMESH1(2),F2),(BMESH1(3),F3),
     > (BMESH1(4),F4),(BMESH1(5),F5),(BMESH1(6),F6),(BMESH1(7),F7),
     >(BMESH1(8),F8),(BMESH1(9),F9)

      SAVE SCAL
      DIMENSION DBDX(3), DBDXI(3)
      SAVE DBDX
      LOGICAL SUMF, SUMFI
      SAVE SUMF
      DIMENSION IQMP(MMAP)
      SAVE IQMP
      SAVE MOD, MOD2
      LOGICAL FIRST

C      DATA IMAP / 1 /       
      DATA DBDX / 3* 0.D0 /
      DATA SUMF / .FALSE. /
      DATA MOD, MOD2 / 0, 0 /
      DATA FIRST / .TRUE. /

C      CALL KSMAP(
C     >           IMAP)

      IF(KUASEX .EQ. 9) THEN
C------- MAP2D, MAP2D-E.  Tracking in symmetryless 2D map
        IF    (KFLD .EQ. MG) THEN
          CALL PAVEL(IRD,A1,R1,Z1,
     >                            B,DB,DDB,KERK)
          CALL NORMB(SCAL,BRI,
     >                        B,DB,DDB)
          BZ=B(1,1)
        ELSEIF(KFLD .EQ. LC) THEN
          CALL PAVEL(IRD,A1,R1,Z1,
     >                            E,DE,DDE,KERK)
          CALL NORMB(SCAL,BRI,
     >                        E,DE,DDE)
        ENDIF
        GOTO 998
      ENDIF

 
      DA=XH(2)-XH(1)
      IAC=INT( (A1-XH(1))/DA+1.5D0 )
      IF(KUASEX .EQ. 8) THEN
C--------  KUASEX = 8:  1-D  map  with  X-cylindrical symmetry (e.g., solenoidal field)
C    1-D  mesh, 5 points grid
          IAC=MAX0(IAC,3)
          IAC=MIN0(IAC,IXMA-2)
          A=A1-XH(IAC)
          GOTO 42
      ENDIF
 
      IF(IRD .EQ. 2) THEN
C-------  2-D/3*3 points or 3-D/3*3*3 points grid
        IAC=MAX0(IAC,2)
        IAC=MIN0(IAC,IXMA-1)
      ELSE
C-------  2-D  5*5 points
        IAC=MAX0(IAC,3)
        IAC=MIN0(IAC,IXMA-2)
      ENDIF
      A=A1-XH(IAC)
 
      DR=YH(2)-YH(1)           
      IRC=INT( (R1-YH(1))/DR+1.5D0 )

C      IF    (IRC.LE.1 .OR. IRC.GE.JYMA) THEN
CC        ... POINT en limite ou HORS CARTE DE Champ
      IF    (IRC.LT.1 .OR. IRC.GT.JYMA) THEN
C        ... POINT HORS CARTE DE Champ
         KERK=1
      ELSE
C        ... Particle is inside field map
         KERK=0
      ENDIF

      IF(IRD .EQ. 2) THEN
C      +++ 2-D  3*3 points or 3-D  3*3*3 points grid
        IRC=MAX0(IRC,2)
        IRC=MIN0(IRC,JYMA-1)
      ELSE
C      +++  2-D grid,  5*5 points
        IRC=MAX0(IRC,3)
        IRC=MIN0(IRC,JYMA-2)
      ENDIF
      R=R1-YH(IRC)
        
      IF(KUASEX .EQ. 7) THEN
C--------  KUASEX = 7: CARTE 3-D   
        IF(IZ.EQ.1) 
     >    CALL ENDJOB('SBR CHAMK, 3D map, give IZ>1 in PARIZ.H',-99)
C        I2=2 introduced to avoid compiler complainig when IZ=1...
        I2 = 2
        DZ=ZH(I2)-ZH(1)
        IZC=INT( (Z1-ZH(1))/DZ+1.5D0 )
        IZC=MAX0(IZC,2)
        IZC=MIN0(IZC,KZMA-1)
        Z=Z1-ZH(IZC)

        IF(MOD.NE.16) THEN
          GOTO 41
        ELSE
          GOTO 51
        ENDIF
      ENDIF

C     .... KUASEX = 1,2,3,4,5,6 : 2-D mid-plane field maps
      IF (IRD.EQ.4) GO TO 21
      IF (IRD.GT.4) GO TO 22
 
      CONTINUE
C       .... ORDRE 2 ,  3*3 points grid

      F1= HC(ID,IAC-1,IRC-1,1,IMAP) * SCAL
      F2= HC(ID,IAC  ,IRC-1,1,IMAP) * SCAL
      F3= HC(ID,IAC+1,IRC-1,1,IMAP) * SCAL
      F4= HC(ID,IAC-1,IRC  ,1,IMAP) * SCAL
      F5= HC(ID,IAC  ,IRC  ,1,IMAP) * SCAL
      F6= HC(ID,IAC+1,IRC  ,1,IMAP) * SCAL
      F7= HC(ID,IAC-1,IRC+1,1,IMAP) * SCAL
      F8= HC(ID,IAC  ,IRC+1,1,IMAP) * SCAL
      F9= HC(ID,IAC+1,IRC+1,1,IMAP) * SCAL

C BB/FM/ 18/08 to be thought of......
      CALL MAPLIM(*999, 9, BMESH1)

      C1=F1+F2+F3+F4+F5+F6+F7+F8+F9
      C2=-F1+F3-F4+F6-F7+F9
      C3=-F1-F2-F3+F7+F8+F9
      C4=F1+F3+F4+F6+F7+F9
      C5=F1+F2+F3+F7+F8+F9
      C6=F1-F3-F7+F9

      BZXXX = 0.D0
      BZXYY = 0.D0
      BZXXX = 0.D0
      BZXYY = 0.D0
      BZXXY = 0.D0
      BZYYY = 0.D0
      BZXXY = 0.D0
      BZYYY = 0.D0
      BZX3Y = 0.D0
      BZXY3 = 0.D0
      BZX3Y = 0.D0
      BZXY3 = 0.D0
      BZX2Y2= 0.D0
      BZY4  = 0.D0
      BZX4  = 0.D0
      BZX4  = 0.D0
      BZX2Y2= 0.D0
      BZY4  = 0.D0

      BZXY=C6/(4.D0*DA*DR)                *BRI
      BZX=C2/(6.D0*DA)                    *BRI
      BZY=C3/(6.D0*DR)                    *BRI
      BZ=(20.D0*C1-12.D0*C4-12.D0*C5)/36.D0   *BRI
      BZXX=(-12.D0*C1+18.D0*C4)/(18.D0*DA*DA) *BRI
      BZYY=(-12.D0*C1+18.D0*C5)/(18.D0*DR*DR) *BRI

c       write(*,*) ' chamk scal, f1-5 ',scal,f1,f2,f3,f4,f5
c     > ,HC(ID,IAC-1,IRC-1,1,IMAP)

      GOTO 30
 
 21   CONTINUE
C       .... ORDRE 4 , GRILLE A 5*5 POINTS
 
C     *** CALCUL DES 15 COEFFS DU POLYNOME DE DEGRE 4 :
C      BZ=A0+A10*X+A11*Y+A20*X2+A21*XY+A22*Y2+...+A42*X2Y2+A43*XY3+A44*Y4
      S1   =0.D0
      SI   =0.D0
      SJ   =0.D0
      SI2  =0.D0
      SIJ  =0.D0
      SJ2  =0.D0
      SI3  =0.D0
      SI2J =0.D0
      SIJ2 =0.D0
      SJ3  =0.D0
      SI4  =0.D0
      SI3J =0.D0
      SI2J2=0.D0
      SIJ3 =0.D0
      SJ4  =0.D0
      DO 1 J=1,5
         JR=3-J
         RJ=DBLE(JR)
         RJ2=RJ*RJ
         RJ3=RJ2*RJ
         IRCJR=IRC+JR
         DO 2 I=1,5
            IA=I-3
            BIAJR= HC(ID,IAC+IA,IRCJR,1,IMAP) * SCAL
            BMESH(I,J) = BIAJR
            AI=DBLE(IA)
            AI2=AI*AI*BIAJR
            AI3=AI2*AI
            AI4=AI3*AI
            AI =AI*BIAJR
            S1   =S1   +         BIAJR
            SI   =SI   +AI
            SJ   =SJ   +RJ      *BIAJR
            SI2  =SI2  +AI2
            SIJ  =SIJ  +AI*RJ
            SJ2  =SJ2  +RJ2     *BIAJR
            SI3  =SI3  +AI3
            SI2J =SI2J +AI2*RJ
            SIJ2 =SIJ2 +AI *RJ2
            SJ3  =SJ3  +RJ3     *BIAJR
            SI4  =SI4  +AI4
            SI3J =SI3J +AI3*RJ
            SI2J2=SI2J2+AI2*RJ2
            SIJ3 =SIJ3 +AI *RJ3
            SJ4  =SJ4  +RJ3*RJ  *BIAJR
2        CONTINUE
1     CONTINUE
 
      CALL MAPLIM(*999, 25, BMESH) 

C     *** CALCUL DE BZ ET SES DERIVEES AU POINT DE MAILLAGE (IAC,IRC)
      DA1  =1.D0/DA
      DA12 =DA1 *DA1
      DA13 =DA12*DA1
      DR1  =1.D0/DR
      DR12 =DR1 *DR1
      DR13 =DR12*DR1
      BZXXX = (-3.4D0 *SI + SI3  )*F1S12
      BZXYY = (-2.D0  *SI + SIJ2 )*F1S70
      BZX   = (.02D0 *SI - 3.4D0*BZXXX*F1S6 - BZXYY) *DA1        *BRI
      BZXXX = BZXXX                                  *DA13       *BRI
      BZXYY = BZXYY                                  *DA1 *DR12  *BRI
      BZXXY = (-2.D0  *SJ + SI2J )*F1S70
      BZYYY = (-3.4D0 *SJ + SJ3  )*F1S12
      BZY   = (.02D0 *SJ - BZXXY - 3.4D0*BZYYY*F1S6)      *DR1   *BRI
      BZXXY = BZXXY                                  *DA12*DR1   *BRI
      BZYYY = BZYYY                                       *DR13  *BRI
      BZX3Y = (-3.4D0 *SIJ+ SI3J )*F1S24
      BZXY3 = (-3.4D0 *SIJ+ SIJ3 )*F1S24
      BZXY  = (.01D0*SIJ- 3.4D0*(BZX3Y + BZXY3)*F1S6)*DA1 *DR1   *BRI
      BZX3Y = BZX3Y                                  *DA13*DR1   *BRI
      BZXY3 = BZXY3                                  *DA1 *DR13  *BRI
      BZX2Y2= (    SI2J2-2.D0*(SI2+SJ2)+ 4.D0 *S1)*F1S49
      BZY4  = (7.D0 *SJ4  -31.D0 *SJ2    +14.4D0 *S1)*F1S12
      BZYY  = ((SJ2-2.D0*S1-310.D0*F1S24*BZY4)*F1S35-BZX2Y2)
      BZX4  = (7.D0 *SI4  -31.D0 *SI2    +14.4D0 *S1)*F1S12
      BZXX  = (SI2-2.D0*S1-155.D0*BZX4*F1S12)*F1S35-BZX2Y2
      BZ =(.01D0*SI2J2 - 1.7D0*(BZXX+BZYY)-13.D0*(BZX4+BZY4)*F1S24 - 
     >2.89D0*BZX2Y2) * BRI
      BZXX  = BZXX   *DA12       *BRI
      BZYY  = BZYY   *DR12       *BRI
      BZX4  = BZX4   *DA13*DA1   *BRI
      BZX2Y2= BZX2Y2 *DA12*DR12  *BRI
      BZY4  = BZY4   *DR13*DR1   *BRI
 
C  CALCUL DE BZ ET SES DERIVEES AU POINT COURANT A1,R1
      A2=A *A
      A3=A2*A
      R2=R *R
      R3=R2*R
      BZ=BZ+BZX*A+BZY*R+0.5D0*BZXX*A2+BZXY*A*R+0.5D0*BZYY*R2
     >+F1S6*BZXXX*A3+.5D0*BZXXY*A2*R 
     > +.5D0*BZXYY*A *R2+F1S6*BZYYY*R3
     >+(0.25D0*BZX4*A+BZX3Y*R)*F1S6*A3+.25D0*BZX2Y2*A2*R2
     >+(BZXY3*A+0.25D0*BZY4*R)*F1S6*R3
      BZX=BZX+  ((BZXX+(BZXXY+.5D0*BZX2Y2*R)*R
     >     +.5D0*(BZXXX+   BZX3Y*R+BZX4 /3.D0*A)*A)*A)
     >       +   (BZXY+.5D0*(BZXYY+  BZXY3/3.D0*R)*R)*R
      BZY=BZY+  ((BZYY+(BZXYY+.5D0*BZX2Y2*A)*A
     >       +.5D0*(BZYYY+   BZXY3*A+BZY4 /3.D0*R)*R)*R)
     >       +   (BZXY+.5D0*(BZXXY+  BZX3Y/3.D0*A)*A)*A
      BZXX=BZXX+(BZXXX+BZX3Y*R+.5D0*BZX4*A)*A
     >         +(BZXXY+.5D0*BZX2Y2*R)*R
      BZXY=BZXY+(BZXXY + .5D0*BZX3Y*A + BZX2Y2*R)*A
     >         +(BZXYY+.5D0*BZXY3*R)*R
      BZYY=BZYY+(BZYYY+BZXY3*A+.5D0*BZY4*R)*R
     >         +(BZXYY+.5D0*BZX2Y2*A)*A
 
CC     ** PASSAGE DE LA FACE MAGNETIQUE (FMAG>.45) ,
CC        SI ON UTILISE CHAMBR
C      IF(LIMIT .EQ. 1) FMAG=BZ*BR/BMAX
 
      IF(KART .NE. 1) THEN
C        ****TRANSFORMATION POLAIRE-CARTESIEN
         R11=1.D0/R1
         R12=R11*R11
         BZX   =   BZX*R11
         BZX4  =(((BZX4 -8.D0*BZXX)*R11+6.D0*BZXXY-3.D0*BZY)*R11+
     >    3.D0*BZYY)*R12
         BZX3Y =(( 6.D0*BZX-3.D0*BZXXX*R11+BZX3Y-8.D0*BZXY)*R11+
     >    3.D0*BZXYY)*R12
         BZX2Y2=(((6.D0*BZXX*R11-4.D0*BZXXY+2.D0*BZY)*R11-
     >    2.D0*BZYY+BZX2Y2)*R11+BZYYY)*R11
         BZXY3 =((6.D0*(BZXY-BZX)*R11-3.D0*BZXYY)*R11+BZXY3)*R11
         BZXXX = ( BZXXX*R11 + 3.D0*BZXY - 2.D0*BZX )*R12
         BZXXY =(( BZXXY - 2.D0*BZXX*R11 - BZY )*R11 + BZYY)*R11
         BZXYY =   BZXYY*R11 + 2.D0*( BZX - BZXY )*R12
         BZXX  = ( BZXX*R11+BZY)*R11
         BZXY  = ( BZXY-BZX)*R11
      ENDIF
 
      GOTO 99
 
 22   CONTINUE
C--------- ORDRE 2 , GRILLE A 5*5 POINTS
 
C     --- CALCUL DES 6 COEFFS DU POLYNOME DE DEGRE 2 :
C       BZ=A0+A10*X+A11*Y+A20*X2+A21*XY+A22*YY

c          write(*,*) ' chamk '
c          write(*,*) ' chamk  SCAL : ',SCAL,DBDX(1),DBDX(2)
      A0 =0D0
      A10=0D0
      A11=0D0
      A20=0D0
      A21=0D0
      A22=0D0
      DO 3 J=1,5
         JR=3-J
         RJ=DBLE(JR)
         RJ2=RJ*RJ
         IRCJR=IRC+JR
         DO 4 I=1,5
            IA=I-3
            AI=DBLE(IA)
            AI2=AI*AI
            BIAJR= HC(ID,IAC+IA,IRCJR,1,IMAP) * SCAL
     >      * (1.D0 + (DBDX(1) + (DBDX(2) + DBDX(3) * RJ) * RJ ) * RJ) 
            A0 =A0 +(27.D0-5.D0*(AI2+RJ2))*BIAJR
            A10=A10+         AI       *BIAJR
            A11=A11+             RJ   *BIAJR
            A20=A20+       ( AI2-2.D0 ) *BIAJR
            A21=A21+         AI *RJ   *BIAJR
            A22=A22+       (-2.D0 +RJ2) *BIAJR
 4       CONTINUE
 3    CONTINUE

      CALL MAPLIM(*999, 25, BMESH)

      FAC=1.D0/70D0
      A0 =A0 *FAC*.4D0
      A10=A10*.02D0
      A11=A11*.02D0
      A20=A20*FAC
      A21=A21*.01D0
      A22=A22*FAC
C
C  CALCUL DE BZ ET SES DERIVEES AU POINT DE MAILLAGE IAC,IRC
      BZXXX = 0.D0
      BZXYY = 0.D0
      BZXXX = 0.D0
      BZXYY = 0.D0
      BZXXY = 0.D0
      BZYYY = 0.D0
      BZXXY = 0.D0
      BZYYY = 0.D0
      BZX3Y = 0.D0
      BZXY3 = 0.D0
      BZX3Y = 0.D0
      BZXY3 = 0.D0
      BZX2Y2= 0.D0
      BZY4  = 0.D0
      BZX4  = 0.D0
      BZX4  = 0.D0
      BZX2Y2= 0.D0
      BZY4  = 0.D0
      BZXX=2.D0*A20/(DA*DA)*BRI
      BZXY=   A21/(DA*DR)  *BRI
      BZYY=2.D0*A22/(DR*DR)*BRI
      BZX=A10    /DA       *BRI
      BZY=A11    /DR       *BRI
      BZ=A0                *BRI
C
 
 30   CONTINUE
C  +++ CALCUL DE BZ ET SES DERIVEES AU POINT COURANT A1,R1
C       dum = BZ
      BZ=BZ+BZX*A+BZY*R+0.5D0*BZXX*A*A+0.5D0*BZYY*R*R+BZXY*A*R
      BZX=BZX+BZXX*A+BZXY*R
      BZY=BZY+BZXY*A+BZYY*R

C         aamm = XH(iac)*180.d0/(4.d0*atan(1.d0))
C        write(88,*) bz*br,dum*br,a,r,iac,aamm,
C     >              HC(ID,IAC  ,IRC  ,1),' chamk,IMAP' 
C        write( *,*) bz*br,dum*br,a1,a,r,iac,aamm,
C     >              HC(ID,IAC  ,IRC  ,1),' chamk,IMAP'


CC     ** PASSAGE DE LA FACE MAGNETIQUE (FMAG>.45) ,
CC        SI ON UTILISE CHAMBR
C      IF(LIMIT .EQ. 1) FMAG=ABS(BZ*BR/BMAX)
 
      IF(KART .NE. 1) THEN
C------- Transformation polaire-cartesien
        BZX=BZX/R1
        BZXXX = (3.D0 * BZXY - 2.D0 * BZX) / (R1 * R1)
        BZXXY = ((-2.D0 * BZXX / R1 - BZY) / R1 + BZYY) / R1
        BZXYY = 2.D0 * (BZX - BZXY) / (R1 * R1)
        BZXX=(BZXX/R1+BZY)/R1
        BZXY=(BZXY-BZX)/R1
      ENDIF
 
      GOTO 99
 
C------------------------------------------------------------------------
 41   CONTINUE
C     .... KUASEX = 7:   3-D  field map (e.g. from TOSCA), 
C     integration ordre 2,   3-D grid with  3*3*3 nodes
 
C   *** Calculation of 6 coeffs of 3 polynomes with degre 2 :
C     BX (Btta)  or  BY (Br)  or  BZ =
C       A000 + A100*X + A010*Y + A001*Z + A200*X2 + A020*Y2 + A002*Z2
C     + A110*XY + A101*XZ + A011*YZ

      DO 415 L = 1,3
        A000(L)=0.D0
        A100(L)=0.D0
        A010(L)=0.D0
        A001(L)=0.D0
        A200(L)=0.D0
        A020(L)=0.D0
        A002(L)=0.D0
        A110(L)=0.D0
        A101(L)=0.D0
        A011(L)=0.D0
        DO 416 J=1,3
          JR=2-J
          DO 416 I=1,3
            IA=I-2
            DO 416 K=1,3
              KZ=K-2
c       if(noel.eq.89)  write(*,*) ' chamk //// ', l,
c     > IAC+IA,IRC+JR,IZC+KZ,IMAP, SCAL
c       if(noel.eq.89)  write(*,*) ' chamk  ', l,
c     > HC(L,IAC+IA,IRC+JR,IZC+KZ,IMAP), SCAL
              BIJK= HC(L,IAC+IA,IRC+JR,IZC+KZ,IMAP) * SCAL
              BMESH3(K,I,J) = BIJK
              A000(L)=A000(L) + 
     >             DBLE(7-3*(IA*IA+JR*JR+KZ*KZ))/3.D0 *BIJK
C              A000(L)=A000(L) +(7.D0/3.D0-IA*IA-JR*JR-KZ*KZ)*BIJK

C                if(IZC.eq.11) then
C                  write(89,fmt='(4i4,2g28.20)') l,j,i,kz,BIJK, A000(L)
C                endif

              A100(L)=A100(L) +        IA        *BIJK
              A010(L)=A010(L) +        JR        *BIJK
              A001(L)=A001(L) +        KZ        *BIJK
              A200(L)=A200(L) +  DBLE(3*IA*IA-2)   *BIJK
              A020(L)=A020(L) +  DBLE(3*JR*JR-2)   *BIJK
              A002(L)=A002(L) +  DBLE(3*KZ*KZ-2)   *BIJK
C              A200(L)=A200(L) +  (3.D0*IA*IA-2.D0)   *BIJK
C              A020(L)=A020(L) +  (3.D0*JR*JR-2.D0)   *BIJK
C              A002(L)=A002(L) +  (3.D0*KZ*KZ-2.D0)   *BIJK
              A110(L)=A110(L) +      IA*JR       *BIJK
              A101(L)=A101(L) +      IA*KZ       *BIJK
              A011(L)=A011(L) +      JR*KZ       *BIJK
 416    CONTINUE

        CALL MAPLIM(*999, 27, BMESH3)


        A000(L)=A000(L)/( 9.D0      )*BRI
        A100(L)=A100(L)/(18.D0*DA   )*BRI
        A010(L)=A010(L)/(18.D0*DR   )*BRI
        A001(L)=A001(L)/(18.D0*DZ   )*BRI
        A200(L)=A200(L)/(18.D0*DA*DA)*BRI
        A020(L)=A020(L)/(18.D0*DR*DR)*BRI
        A002(L)=A002(L)/(18.D0*DZ*DZ)*BRI
        A110(L)=A110(L)/(12.D0*DA*DR)*BRI
        A101(L)=A101(L)/(12.D0*DA*DZ)*BRI
        A011(L)=A011(L)/(12.D0*DR*DZ)*BRI
 415  CONTINUE
C
C  CALCUL BZ ET SES DERIVEES AU POINT COURANT A1,R1,Z1
C
C     *** COMPOSANTES BX, BY, BZ DU Champ

      DO 417 L = 1,3
        B(1,L)=A000(L)     
     >       + A100(L)*A   + A010(L)*R   + A001(L)*Z
     >       + A200(L)*A*A + A020(L)*R*R + A002(L)*Z*Z
     >       + A110(L)*A*R + A101(L)*A*Z + A011(L)*R*Z
        DB(1,L)  = A100(L) + 2.D0*A200(L)*A + A110(L)*R + A101(L)*Z
        DB(2,L)  = A010(L) + 2.D0*A020(L)*R + A110(L)*A + A011(L)*Z
        DB(3,L)  = A001(L) + 2.D0*A002(L)*Z + A101(L)*A + A011(L)*R
        DDB(1,1,L) = 2.D0*A200(L)
        DDB(1,2,L) = A110(L)
        DDB(2,1,L) = A110(L)
        DDB(2,2,L) = 2.D0*A020(L)
        DDB(1,3,L) = A101(L)
        DDB(3,1,L) = A101(L)
        DDB(2,3,L) = A011(L)
        DDB(3,2,L) = A011(L)
        DDB(3,3,L) = 2.D0*A002(L)
 417  CONTINUE
      BZ = B(1,3)

C-------------- Pour éventuel tests si défaut de plan médian 
C          if(izc.eq.21) then   ! KEK FFAG
C  TEST RACCAM
C          if(izc.eq.11) then   ! RACCAM
C            write(89,*) r1,z1,b(1,1), b(1,2), ' sbr chamk'
C            b(1,1)=0.d0
C            b(1,2)=0.d0
C            DB(1,1)  = 0.d0
C            DB(2,1)  = 0.d0
C            dDB(1,1,1)  = 0.d0
C            dDB(1,1,2)  = 0.d0
C            dDB(1,2,1)  = 0.d0
C            dDB(1,2,2)  = 0.d0
C            dDB(2,1,1)  = 0.d0
C            dDB(2,1,2)  = 0.d0
C          endif
C-----------------------------------------------


C For calculation of FF coefficients in DFD ffag triplet using mathematica
C         if(a1.le.0.d0) write(88,fmt='(1p,3g14.6)') a1,r1,bz*br

CC------- PASSAGE DE LA FACE MAGNETIQUE (FMAG>.45) ,
CC        SI ON UTILISE CHAMBR
C      IF(LIMIT .EQ. 1) FMAG=ABS(BZ*BR/BMAX)
 
      IF(KART .NE. 1) THEN
C--------- Transformation from cylindrical to cartesian coordinates
        R11=1.D0/R1
        R12=R11*R11
        DO 418 L = 1,3
C         dB_*/dtta -> /rdB_*/dx 
          DB(1,L)  = DB(1,L)*R11
C         d2B_*/dtta2 -> d2B_*/dx2
          DDB(1,1,L) = ( DDB(1,1,L)*R11 + DB(2,L) ) * R11
C         d2B_*/dtta.dr -> d2B_*/dxdy
          DDB(1,2,L) = ( DDB(1,2,L) - DB(1,L) ) * R11
          DDB(2,1,L) = DDB(1,2,L) 
C         d2B_*/dtta.dz -> d2B_*/dxdz
          DDB(1,3,L) = DDB(1,3,L) * R11
          DDB(3,1,L) = DDB(1,3,L)
C--------- Call to DBDXYZ is redundant because all derivatives 
C          have been filled above, yet it links the derivatives by Maxwell 
          IDB=2
C          CALL DBDXYZ(IDB,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
 418    CONTINUE
      ENDIF
 
      GOTO 998
 
 
C------------------------------------------------------------------------
 51   CONTINUE
C     .... KUASEX = 7 and MOD=16 : same as MOD.ne.16, apart from the following : 
C     Field contributions from all maps are summed right here, whereas otherwise field 
C     maps are summed up into a single new one when they are read (managed in toscac.f). 

      DO 515 L = 1,3
        A000(L)=0.D0
        A100(L)=0.D0
        A010(L)=0.D0
        A001(L)=0.D0
        A200(L)=0.D0
        A020(L)=0.D0
        A002(L)=0.D0
        A110(L)=0.D0
        A101(L)=0.D0
        A011(L)=0.D0
        DO J=1,3
          JR=2-J
          DO I=1,3
            IA=I-2
            DO K=1,3
              KZ=K-2
              BIJK= HC(L,IAC+IA,IRC+JR,IZC+KZ,JMAP) * SCAL
              BMESH3(K,I,J) = BIJK
              A000(L)=A000(L) + 
     >             DBLE(7-3*(IA*IA+JR*JR+KZ*KZ))/3.D0 *BIJK
              A100(L)=A100(L) +        IA        *BIJK
              A010(L)=A010(L) +        JR        *BIJK
              A001(L)=A001(L) +        KZ        *BIJK
              A200(L)=A200(L) +  DBLE(3*IA*IA-2)   *BIJK
              A020(L)=A020(L) +  DBLE(3*JR*JR-2)   *BIJK
              A002(L)=A002(L) +  DBLE(3*KZ*KZ-2)   *BIJK
              A110(L)=A110(L) +      IA*JR       *BIJK
              A101(L)=A101(L) +      IA*KZ       *BIJK
              A011(L)=A011(L) +      JR*KZ       *BIJK
            ENDDO
          ENDDO
        ENDDO

        CALL MAPLIM(*999, 27, BMESH3)

        A000(L)=A000(L)/( 9.D0      )*BRI
        A100(L)=A100(L)/(18.D0*DA   )*BRI
        A010(L)=A010(L)/(18.D0*DR   )*BRI
        A001(L)=A001(L)/(18.D0*DZ   )*BRI
        A200(L)=A200(L)/(18.D0*DA*DA)*BRI
        A020(L)=A020(L)/(18.D0*DR*DR)*BRI
        A002(L)=A002(L)/(18.D0*DZ*DZ)*BRI
        A110(L)=A110(L)/(12.D0*DA*DR)*BRI
        A101(L)=A101(L)/(12.D0*DA*DZ)*BRI
        A011(L)=A011(L)/(12.D0*DR*DZ)*BRI
 515  CONTINUE
C
C  CALCUL BZ ET SES DERIVEES AU POINT COURANT A1,R1,Z1
C
C     *** COMPOSANTES BX, BY, BZ DU Champ

      DO 517 L = 1,3
        B(1,L)=A000(L)     
     >       + A100(L)*A   + A010(L)*R   + A001(L)*Z
     >       + A200(L)*A*A + A020(L)*R*R + A002(L)*Z*Z
     >       + A110(L)*A*R + A101(L)*A*Z + A011(L)*R*Z
        DB(1,L)  = A100(L) + 2.D0*A200(L)*A + A110(L)*R + A101(L)*Z
        DB(2,L)  = A010(L) + 2.D0*A020(L)*R + A110(L)*A + A011(L)*Z
        DB(3,L)  = A001(L) + 2.D0*A002(L)*Z + A101(L)*A + A011(L)*R
        DDB(1,1,L) = 2.D0*A200(L)
        DDB(1,2,L) = A110(L)
        DDB(2,1,L) = A110(L)
        DDB(2,2,L) = 2.D0*A020(L)
        DDB(1,3,L) = A101(L)
        DDB(3,1,L) = A101(L)
        DDB(2,3,L) = A011(L)
        DDB(3,2,L) = A011(L)
        DDB(3,3,L) = 2.D0*A002(L)
 517  CONTINUE
      BZ = B(1,3)

C-------------- Pour éventuel tests si défaut de plan médian 
C          if(izc.eq.21) then   ! KEK FFAG
C  TEST RACCAM
C          if(izc.eq.11) then   ! RACCAM
C            write(89,*) r1,z1,b(1,1), b(1,2), ' sbr chamk'
C            b(1,1)=0.d0
C            b(1,2)=0.d0
C            DB(1,1)  = 0.d0
C            DB(2,1)  = 0.d0
C            dDB(1,1,1)  = 0.d0
C            dDB(1,1,2)  = 0.d0
C            dDB(1,2,1)  = 0.d0
C            dDB(1,2,2)  = 0.d0
C            dDB(2,1,1)  = 0.d0
C            dDB(2,1,2)  = 0.d0
C          endif
C-----------------------------------------------


C For calculation of FF coefficients in DFD ffag triplet using mathematica
C         if(a1.le.0.d0) write(88,fmt='(1p,3g14.6)') a1,r1,bz*br

CC------- PASSAGE DE LA FACE MAGNETIQUE (FMAG>.45) ,
CC        SI ON UTILISE CHAMBR
C      IF(LIMIT .EQ. 1) FMAG=ABS(BZ*BR/BMAX)
 
      IF(KART .NE. 1) THEN
C--------- Transformation from cylindrical to cartesian coordinates
        R11=1.D0/R1
        R12=R11*R11
        DO L = 1,3
C         dB_*/dtta -> /rdB_*/dx 
          DB(1,L)  = DB(1,L)*R11
C         d2B_*/dtta2 -> d2B_*/dx2
          DDB(1,1,L) = ( DDB(1,1,L)*R11 + DB(2,L) ) * R11
C         d2B_*/dtta.dr -> d2B_*/dxdy
          DDB(1,2,L) = ( DDB(1,2,L) - DB(1,L) ) * R11
          DDB(2,1,L) = DDB(1,2,L) 
C         d2B_*/dtta.dz -> d2B_*/dxdz
          DDB(1,3,L) = DDB(1,3,L) * R11
          DDB(3,1,L) = DDB(1,3,L)
        ENDDO
        IDB=2
      ENDIF
 
      GOTO 998
 
 
 42   CONTINUE
C------- KUASEX = 8: 1-D map with  x-revolution symmetry, grid 1-D with 5 points
        CALL KRTAX(A1,BX,SCAL)
        R2=R1*R1+Z1*Z1
        R=SQRT(R2)
        CALL BAXBXR(BX,R,R2,BC,DBC,DDBC)
        IF    (KFLD .EQ. MG) THEN
C--------- magnetique
          CALL BXRXYZ(BC,DBC,DDBC,R1,Z1,R,2,B,DB,DDB)
          CALL DBDXYZ(2,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
        ELSEIF(KFLD .EQ. LC) THEN
C--------- electrique
          CALL BXRXYZ(BC,DBC,DDBC,R1,Z1,R,2,E,DE,DDE)
          CALL DBDXYZ(2,DE,DDE,D3EX,D3EY,D3EZ,D4EX,D4EY,D4EZ)
        ENDIF
      GOTO 998
 
 
 99   CONTINUE

C Due to reminicences from ancient times
      BZ0(1,1)=BZ
      BZ0(2,1)=BZX
      BZ0(1,2)=BZY
      BZ0(3,1)=BZXX
      BZ0(2,2)=BZXY
      BZ0(1,3)=BZYY
      BZ0(4,1)=BZXXX
      BZ0(3,2)=BZXXY
      BZ0(2,3)=BZXYY
      BZ0(1,4)=BZYYY
      BZ0(5,1)=BZX4
      BZ0(4,2)=BZX3Y
      BZ0(3,3)=BZX2Y2
      BZ0(2,4)=BZXY3
      BZ0(1,5)=BZY4
 
      IDZ=2
      IF(KART .NE. 1) IDZ=3
      IF(IRD .EQ. 4) IDZ=4
      CALL SYMMED(Z1,IDZ,BZ0,
     >                       B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY)
      IDB=2
      IF(KART .NE. 1) IDB=3
      IF(IRD .EQ. 4) IDB=4
      CALL DBDXYZ(IDB,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
 
 998  CONTINUE

      RETURN

 999  RETURN 1

      ENTRY CHAMK2(SCALI)
      SCAL = SCALI
      RETURN

      ENTRY CHAMK4(DBDXI,IND)
      DO I = 1, IND
        DBDX(I) = DBDXI(I)
      ENDDO
      RETURN

      ENTRY CHAMK6(SUMFI,MODI,MOD2I)
C MOD2 is the number of field maps under this TOSCA[MOD.MOD2=16.MOD2] set
C IMAP is the current number for the first map in this TOSCA[MOD.MOD2=16.MOD2] set
      SUMF = SUMFI
      MOD = MODI
      MOD2 = MOD2I
      CALL KSMAP(
     >           IMAP)
      IF(IMAP+MOD2-1 .GT.MMAP) CALL ENDJOB('Pgm chamk. Too many field' 
     >//' maps under this TOSCA keyword. Max allowed is ',MMAP-IMAP+1)
      DO I = 1, MOD2
        IQMP(I) = IMAP+I-1
      ENDDO
      FIRST = .TRUE.
      RETURN

      END
