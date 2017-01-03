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
C  -------
      SUBROUTINE INTPLF(R1,A,R,DA,DR,FTAB,IRD, 
     >                                           BZ0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     -----------------------------
C     Interpolation of field values
C     -----------------------------
      DIMENSION FTAB(5,5),BZ0(5,5)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH
 
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      INCLUDE "C.MARK.H"     ! COMMON/MARK/ KART,KALC,KERK,KUASEX
C      LOGICAL ZSYM
      INCLUDE "C.TYPFLD.H"     ! COMMON/TYPFLD/ KFLD,MG,LC,ML,ZSYM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
 
      PARAMETER (F1S6=1.D0/6.D0, F1S12=1.D0/12.D0, F1S24=1.D0/24.D0,
     >F1S35=1.D0/35.D0, F1S49=1.D0/49.D0, F1S70=1.D0/70.D0)
 
      IF(IRD.EQ.2) THEN
C       .... ORDRE 2 , uses a 3*3 points grid

        F1=FTAB(1,3)      !H(ID,IAC-1,IRC-1,1)    FTAB(TTA,RO) : 
        F2=FTAB(2,3)      !H(ID,IAC  ,IRC-1,1)      11 21 31    11 : TTAMIN-RMAX
        F3=FTAB(3,3)      !H(ID,IAC+1,IRC-1,1)      12 22 32
        F4=FTAB(1,2)      !H(ID,IAC-1,IRC  ,1)      13 23 33          33 : TTAMAX-RMIN
        F5=FTAB(2,2)      !H(ID,IAC  ,IRC  ,1)
        F6=FTAB(3,2)      !H(ID,IAC+1,IRC  ,1)
        F7=FTAB(1,1)      !H(ID,IAC-1,IRC+1,1)
        F8=FTAB(2,1)      !H(ID,IAC  ,IRC+1,1)
        F9=FTAB(3,1)      !H(ID,IAC+1,IRC+1,1)

C        write(*,*) ' intplf F1-F4 : ',f1,f2,f3,f4
C        write(*,*) ' intplf F5-F9 : ',f5,f6,f7,f7,f9

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

        IF(A.NE.0.D0 .OR. R.NE.0.D0) THEN
          BZ=BZ+BZX*A+BZY*R+0.5D0*BZXX*A*A+0.5D0*BZYY*R*R+BZXY*A*R
          BZX=BZX+BZXX*A+BZXY*R
          BZY=BZY+BZXY*A+BZYY*R
        ENDIF

      ELSEIF(IRD.EQ.4) THEN
C       .... ORDRE 4 , USES A  5*5 POINTS GRID. 
C      CALCUL DES 15 COEFFS DU POLYNOME DE DEGRE 4 :
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
           DO 2 I=1,5
              IA=I-3
              BIAJR=FTAB(I,J)               ! 11  21 31 41 51,  12  22 32...
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
2          CONTINUE
1       CONTINUE
 
C     *** CALCUL DE BZ ET SES DERIVEES AU POINT COURANT
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
     >  2.89D0*BZX2Y2) * BRI
        BZXX  = BZXX   *DA12       *BRI
        BZYY  = BZYY   *DR12       *BRI
        BZX4  = BZX4   *DA13*DA1   *BRI
        BZX2Y2= BZX2Y2 *DA12*DR12  *BRI
        BZY4  = BZY4   *DR13*DR1   *BRI
 
        IF(A.NE.0.D0 .OR. R.NE.0.D0) THEN
          A2=A *A
          A3=A2*A
          R2=R *R
          R3=R2*R
          BZ=BZ+BZX*A+BZY*R+0.5D0*BZXX*A2+BZXY*A*R+0.5D0*BZYY*R2
     >    +F1S6*BZXXX*A3+.5D0*BZXXY*A2*R 
     >    +.5D0*BZXYY*A *R2+F1S6*BZYYY*R3
     >    +(0.25D0*BZX4*A+BZX3Y*R)*F1S6*A3+.25D0*BZX2Y2*A2*R2
     >    +(BZXY3*A+0.25D0*BZY4*R)*F1S6*R3
          BZX=BZX+  ((BZXX+(BZXXY+.5D0*BZX2Y2*R)*R
     >     +.5D0*(BZXXX+   BZX3Y*R+BZX4 /3.D0*A)*A)*A)
     >         +   (BZXY+.5D0*(BZXYY+  BZXY3/3.D0*R)*R)*R
          BZY=BZY+  ((BZYY+(BZXYY+.5D0*BZX2Y2*A)*A
     >         +.5D0*(BZYYY+   BZXY3*A+BZY4 /3.D0*R)*R)*R)
     >         +   (BZXY+.5D0*(BZXXY+  BZX3Y/3.D0*A)*A)*A
          BZXX=BZXX+(BZXXX+BZX3Y*R+.5D0*BZX4*A)*A
     >           +(BZXXY+.5D0*BZX2Y2*R)*R
          BZXY=BZXY+(BZXXY + .5D0*BZX3Y*A + BZX2Y2*R)*A
     >           +(BZXYY+.5D0*BZXY3*R)*R
          BZYY=BZYY+(BZYYY+BZXY3*A+.5D0*BZY4*R)*R
     >           +(BZXYY+.5D0*BZX2Y2*A)*A
        ENDIF

      ELSEIF(IRD.EQ. 25) THEN
C--------- ORDRE 2 , GRILLE A 5*5 POINTS
 
C     --- CALCUL DES 6 COEFFS DU POLYNOME DE DEGRE 2 :
C       BZ=A0+A10*X+A11*Y+A20*X2+A21*XY+A22*YY
        A0 =0.D0
        A10=0.D0
        A11=0.D0
        A20=0.D0
        A21=0.D0
        A22=0.D0
        DO 3 J=1,5
           JR=3-J
           RJ=DBLE(JR)
           RJ2=RJ*RJ
C           IRCJR=IRC+JR
           DO 4 I=1,5
              IA=I-3
              AI=DBLE(IA)
              AI2=AI*AI
C              BIAJR=H(ID,IAC+IA,IRCJR,1)
              BIAJR=FTAB(I,J) 
              A0 =A0 +(27.D0-5.D0*(AI2+RJ2))*BIAJR
              A10=A10+         AI       *BIAJR
              A11=A11+             RJ   *BIAJR
              A20=A20+       ( AI2-2.D0 ) *BIAJR
              A21=A21+         AI *RJ   *BIAJR
              A22=A22+       (-2.D0 +RJ2) *BIAJR
 4         CONTINUE
 3      CONTINUE

        FAC=1.D0/70.D0
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
 
        IF(A.NE.0.D0 .OR. R.NE.0.D0) THEN
          BZ=BZ+BZX*A+BZY*R+0.5D0*BZXX*A*A+0.5D0*BZYY*R*R+BZXY*A*R
          BZX=BZX+BZXX*A+BZXY*R
          BZY=BZY+BZXY*A+BZYY*R
        ENDIF

      ENDIF

CC     ** PASSAGE DE LA FACE MAGNETIQUE (FMAG>.45) ,
CC        SI ON UTILISE CHAMBR
C      IF(LIMIT .EQ. 1) FMAG=ABS(BZ*BR/BMAX)
 
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
Cccc BZXYY will yield BZ0(2,3) and then B122 and D3BX(3,2,2) = B122 in symmed
         BZXYY =   BZXYY*R11 + 2.D0*( BZX - BZXY )*R12
         BZXX  = ( BZXX*R11+BZY)*R11
         BZXY  = ( BZXY-BZX)*R11
      ENDIF

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
 
      RETURN
      END
