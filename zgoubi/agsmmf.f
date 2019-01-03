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
      SUBROUTINE AGSMMF(ID,X,Y0,Z0,BM,DLEM,DLSM,DE,DS,RT
     >,XE,XS,CE,CS
     >,B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ,BT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(MPOL=10)
      PARAMETER(MCOEF=6)
      DIMENSION BM(MPOL),DLEM(MPOL),DLSM(MPOL)
     >,DE(MPOL,MCOEF),DS(MPOL,MCOEF),RT(MPOL)
      DIMENSION CE(MCOEF),CS(MCOEF)
      DIMENSION B(5,3),DB(3,3),DDB(3,3,3)
      DIMENSION D3BX(3,3,3), D3BY(3,3,3), D3BZ(3,3,3)
      DIMENSION D4BX(3,3,3,3) ,D4BY(3,3,3,3) ,D4BZ(3,3,3,3)
      DIMENSION BT(5, *)

      INCLUDE "C.AIM.H"     ! COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH

      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN
      INCLUDE "C.INTEG.H"     ! COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,HDPRF,DP,QBR,BRI

      LOGICAL CHFE, CHFS, CHU
      PARAMETER (Q7=3.038194444D-4)

      XLS=XLIM-XS
      IOP=1

      IF( BM(1) .EQ. ZERO) GOTO 82

C----- DIPOLE
      IP=1
      GN = BM(IP)*BRI
      DLE = DLEM(IP)
      DLS = DLSM(IP)
      Y=Y0
      Z=Z0
      IF(RT(IP) .NE. ZERO) CALL ROTX(RT(IP),Y,Z)

      IF(DLE+DLS .EQ. ZERO) THEN
C------- DIPOLE CRENEAU

C        IF ( X.GE.XE .AND. X.LE.XS ) THEN
C        IF ( X.GE.XE .AND. X.LT.XS ) THEN
        IF ( X.GE.XE .AND. X.LE.XS ) THEN
C         *** A L'INTERIEUR DE L'AIMANT

          B(1,3)  = GN

        ENDIF

      ELSEIF(DLE+DLS .NE. ZERO) THEN
C       *** DIPOLE AVEC ChampS DE FUITE D'ENTREE ET/OU DE SORTIE

C       DISTANCE A LA FACE D'ENTREE :
        SE=XE-X
C       DISTANCE A LA FACE DE SORTIE :
        SS=X-XS

        ABSE=SE
        IF(ABSE.LT.ZERO) ABSE=-ABSE
        ABSS=SS
        IF(ABSS.LT.ZERO) ABSS=-ABSS

        CHFE=ABSE.LE.XE .AND. DLE .NE. ZERO
        CHU= ABSE.GT.XE .AND. ABSS.GT.XLS
        CHFS=ABSS.LE.XLS .AND. DLS .NE. ZERO

        IF(CHU) THEN
C         *** POINT DANS ZONE DE Champ Uniforme ( G(X)=0 )

            B(1,3)  = GN

        ELSE
C       *** POINT DANS ZONE DE Champ DE FUITE D'ENTRE ET/OU DE SORTIE

          IF(.NOT. CHFE) THEN
C           *** HORS D'EFFET DU Champ DE FUITE D'ENTREE
C             OU GRADIENT CRENEAU

            GE  =UN
            DGE =ZERO
            D2GE=ZERO
            D3GE=ZERO
            D4GE=ZERO
            D5GE=ZERO
            D6GE=ZERO

          ELSE
C           *** EFFET DU Champ DE FUITE D'ENTREE

            SE=SE/DLE
            CALL DRVG(4,CE,SE,GE,DGE,D2GE,D3GE,D4GE,D5GE,D6GE)
C            CALL DRVG(6,CE,SE,GE,DGE,D2GE,D3GE,D4GE,D5GE,D6GE)

          ENDIF

          IF(.NOT. CHFS) THEN
C           *** HORS D'EFFET DU Champ DE FUITE DE SORTIE
C             OU GRADIENT CRENEAU

            GS  =UN
            DGS =ZERO
            D2GS=ZERO
            D3GS=ZERO
            D4GS=ZERO
            D5GS=ZERO
            D6GS=ZERO

          ELSE
C           *** EFFET DU Champ DE FUITE DE SORTIE

            SS=SS/DLS
            CALL DRVG(4,CS,SS,GS,DGS,D2GS,D3GS,D4GS,D5GS,D6GS)
C            CALL DRVG(6,CS,SS,GS,DGS,D2GS,D3GS,D4GS,D5GS,D6GS)

          ENDIF

          G= GE+GS-UN

          IF(G.LT.ZERO) CALL ENDJOB(
     >      'SBR CHAMC : gradient is wrong,  G (=GE+GS-1) < ',0)

C         *** CACULE LE DE Champ DE FUITE RESULTANT G=GE+GS-1, ET
C             SES DERIVEES
          G   =  G*GN
          DG  = (DGE  * DE(IP,1) + DGS  * DS(IP,1))  *BRI
          D2G = (D2GE * DE(IP,2) + D2GS * DS(IP,2))  *BRI
          D3G = (D3GE * DE(IP,3) + D3GS * DS(IP,3))  *BRI
          D4G = (D4GE * DE(IP,4) + D4GS * DS(IP,4))  *BRI
          D5G = (D5GE * DE(IP,5) + D5GS * DS(IP,5))  *BRI
          D6G = (D6GE * DE(IP,6) + D6GS * DS(IP,6))  *BRI

          Y2 = Y  * Y
          Y4 = Y2 * Y2
          Y6 = Y4 * Y2

          Z2 = Z  * Z
          Z4 = Z2 * Z2
          Z6 = Z4 * Z2

          YZ   = Y  * Z

          Y2Z2 = Y2 * Z2

          YPZ= Y2+Z2

          Y4YZZ4=Y4+6.D0*Y2Z2+5.D0*Z4

          Y2P3Z2=Y2+3.D0*Z2

C--------- Components Bx, By, Bz of field
          B(1,1)= Z*(DG + (- D3G/8.D0 + D5G/192.D0*YPZ)*YPZ)
          B(1,2)=  YZ*(-D2G/4.D0 + (D4G/48.D0  - D6G/1536.D0*YPZ )*YPZ)
          B(1,3)=   G - D2G/8.D0*Y2P3Z2 + D4G/192.D0*Y4YZZ4
     >     - D6G/9216.D0*(Y6+9.D0*Y4*Z2+15.D0*Y2*Z4+7.D0*Z6)

c     >  D6G/9216.D0,Y2P3Z2,Y4YZZ4
C         ... dBx/dX
          DB(1,1) = Z*(D2G+( -D4G/8.D0+D6G/192.D0*YPZ )*YPZ)
C         ... dBx/dY = dBy/dX
          DB(2,1) = YZ*(-.25D0*D3G + D5G/48.D0*YPZ)
C         ... dBy/dY
          DB(2,2) = Z*(-.25D0*D2G + D4G/48.D0*(3.D0*Y2+Z2)
     >     - D6G/1536.D0*(5.D0*Y4+6.D0*Y2Z2+Z4))
C     >     - D6G/256.D0*(5.D0*Y4+6.D0*Y2Z2+Z4))
C         .. dBx/dZ = dBz/dX
          DB(3,1) = DG - D3G/8.D0*Y2P3Z2 + D5G/192.D0*Y4YZZ4
C         .. dBy/dZ = dBZ/dY
          DB(3,2) = Y*(-.25D0*D2G + D4G/48.D0*Y2P3Z2-D6G/1536.D0*Y4YZZ4)

C         ... d2Bx/dX2
          DDB(1,1,1) = Z*(D3G - D5G/8.D0*YPZ)
C         ... d2Bx/dXdY = d2By/dX2
          DDB(2,1,1)= YZ*(-.25D0*D4G + D6G/48.D0*YPZ)
C         ... d2Bx/dY2
          DDB(2,2,1) = Z*(-.25D0*D3G + D5G/48.D0*(3.D0*Y2+Z2))
C         ... d2By/dY2
          DDB(2,2,2) = YZ*(D4G/8.D0 - D6G/384.D0*(5.D0*Y2+3.D0*Z2))
C          DDB(2,2,2) = YZ*(D4G/8.D0 - D6G/64.D0*(5.D0*Y2+3.D0*Z2))
C         .. d2Bx/dXdZ = d2Bz/dX2
          DDB(3,1,1) = D2G-D4G/8.D0*Y2P3Z2+D6G/192.D0*Y4YZZ4
C         .. d2By/dXdZ = d2Bz/dXdY = d2Bx/dYdZ
          DDB(3,2,1) = Y*(-.25D0*D3G + D5G/48.D0*Y2P3Z2)
C         .. d2By/dYdZ = d2Bz/dY2
          DDB(3,2,2) = -.25D0*D2G+D4G/16.D0*YPZ-D6G/1536.D0*
     >     (5.D0*Y4+18.D0*Y2Z2+5.D0*Z4)

          IF(ID .GE. 3) THEN
C--------- THIRD AND FOURTH ORDER NOT INSTALLED
          ENDIF
C         ... ENDIF ID>=3

        ENDIF
C         endif CHU

      ENDIF
C     **  endif test DLE+DLS

      IF(RT(IP) .NE. ZERO)
     >  CALL XROTB(RT(IP),B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
C-----  END DIP

      CALL ADPOL(ID,IOP,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)

 82   CONTINUE

      IF( BM(2) .EQ. ZERO) GOTO 83

C--------------------------------------------------------------------
C----- QUADRUPOLE
      IP=2
      GN = BM(IP)*BRI
      DLE = DLEM(IP)
      DLS = DLSM(IP)
      Y=Y0
      Z=Z0

      IF(RT(IP) .NE. ZERO) CALL ROTX(RT(IP),Y,Z)

      IF(DLE+DLS .EQ. ZERO) THEN
C       *** Q-POLE CRENEAU : PAS DE Champ DE FUITE

C                          X.LT.XS rather than X.LE.XS ensures that
C                          the pole contributes zero beyond exit EFB in mixed Sharp Edge + FF
C        IF ( X.GE.XE .AND. X.LT.XS ) THEN
        IF ( X.GE.XE .AND. X.LE.XS ) THEN
C         *** A L'INTERIEUR DE L'AIMANT

          B(1,2)  = GN * Z
          B(1,3)  = GN * Y

          DB(3,2) = GN

        ENDIF

      ELSEIF(DLE+DLS .NE. ZERO) THEN
C     *** QUADRUPOLE AVEC ChampS DE FUITE D'ENTREE ET/OU DE SORTIE

C     DISTANCE A LA FACE D'ENTREE :
      SE=XE-X
C     DISTANCE A LA FACE DE SORTIE :
      SS=X-XS

      ABSE=SE
      IF(ABSE.LT.ZERO) ABSE=-ABSE
      ABSS=SS
      IF(ABSS.LT.ZERO) ABSS=-ABSS

      CHFE=ABSE.LE.XE .AND. DLE .NE. ZERO
      CHU= ABSE.GT.XE .AND. ABSS.GT.XLS
      CHFS=ABSS.LE.XLS .AND. DLS .NE. ZERO

      IF(CHU) THEN
C       *** POINT DANS ZONE DE Champ Uniforme ( G(X)=0 )

          B(1,2)  = GN * Z
          B(1,3)  = GN * Y

          DB(3,2) = GN
      ELSE
C     *** POINT DANS ZONE DE Champ DE FUITE D'ENTRE ET/OU DE SORTIE

        IF(.NOT. CHFE) THEN
C         *** HORS D'EFFET DU Champ DE FUITE D'ENTREE
C             OU GRADIENT CRENEAU

          GE  =UN
          DGE =ZERO
          D2GE=ZERO
          D3GE=ZERO
          D4GE=ZERO
          D5GE=ZERO
          D6GE=ZERO

        ELSE
C         *** EFFET DU Champ DE FUITE D'ENTREE

          SE=SE/DLE
          CALL DRVG(4,CE,SE,GE,DGE,D2GE,D3GE,D4GE,D5GE,D6GE)
C          CALL DRVG(6,CE,SE,GE,DGE,D2GE,D3GE,D4GE,D5GE,D6GE)

        ENDIF

        IF(.NOT. CHFS) THEN
C         *** HORS D'EFFET DU Champ DE FUITE DE SORTIE
C             OU GRADIENT CRENEAU

          GS  =UN
          DGS =ZERO
          D2GS=ZERO
          D3GS=ZERO
          D4GS=ZERO
          D5GS=ZERO
          D6GS=ZERO

        ELSE
C         *** EFFET DU Champ DE FUITE DE SORTIE

          SS=SS/DLS
          CALL DRVG(4,CS,SS,GS,DGS,D2GS,D3GS,D4GS,D5GS,D6GS)
C          CALL DRVG(6,CS,SS,GS,DGS,D2GS,D3GS,D4GS,D5GS,D6GS)

        ENDIF

        G= GE+GS-UN

        IF(G.LT.ZERO) CALL ENDJOB(
     >    'SBR CHAMC : gradient  is  wrong  G (=GE+GS-1) < ',0)

C       *** CACULE LE DE Champ DE FUITE RESULTANT G=GE+GS-1, ET
C           SES DERIVEES
      G   =  G*GN
      DG  = (DGE  * DE(IP,1) + DGS  * DS(IP,1))  *BRI
      D2G = (D2GE * DE(IP,2) + D2GS * DS(IP,2))  *BRI
      D3G = (D3GE * DE(IP,3) + D3GS * DS(IP,3))  *BRI
      D4G = (D4GE * DE(IP,4) + D4GS * DS(IP,4))  *BRI
      D5G = (D5GE * DE(IP,5) + D5GS * DS(IP,5))  *BRI
      D6G = (D6GE * DE(IP,6) + D6GS * DS(IP,6))  *BRI

      Y2 = Y  * Y
      Y3 = Y2 * Y
      Y4 = Y3 * Y
      Y5 = Y4 * Y
      Y6 = Y5 * Y

      Z2 = Z  * Z
      Z3 = Z2 * Z
      Z4 = Z3 * Z
      Z5 = Z4 * Z
      Z6 = Z5 * Z

      YZ   = Y  * Z
      YZ2  = Y  * Z2
      YZ3  = Y  * Z3
      YZ4  = Y  * Z4
      YZ5  = Y  * Z5

      Y2Z2 = Y2 * Z2
      Y2Z3 = Y2 * Z3

      Y3Z  = Y3 * Z
      Y3Z2 = Y3 * Z2
      Y3Z3 = Y3 * Z3

C     *** Bx, By, Bz field components
      B(1,1)=
     > YZ*(DG -D3G/12.D0*(Y2+ Z2)+D5G/384.D0*(Y4+2.D0*Y2Z2+Z4))
      B(1,2)=
     > Z*(G - D2G/12.D0*(Z2+3.D0*Y2) + D4G/384.D0*(Z4+6.D0*Y2Z2+
     > 5.D0*Y4) - D6G/23040.D0*(Z6+9.D0*Y2*Z4+15.D0*Z2*Y4+7.D0*Y6))
      B(1,3)=
     > Y*(G - D2G/12.D0*(Y2+3.D0*Z2) + D4G/384.D0*(Y4+6.D0*Y2Z2+
     > 5.D0*Z4) - D6G/23040.D0*(Y6+9.D0*Y4*Z2+15.D0*Y2*Z4+7.D0*Z6))

C     *** DERIVEES DANS LE PLAN MEDIAN
      BZX   =  DG*Y - D3G*Y3/12.D0  +   D5G*Y5/384.D0
      BZY   = (  G   - D2G*Y2*.25D0 +   D4G*Y4/76.8D0   -  D6G*Y6*Q7 )
      BZXX  = D2G*Y - D4G*Y3/12.D0  +   D6G*Y5/384.D0
      BZXY  =  DG   - D3G*Y2*.25D0 +   D5G*Y4/76.8D0
      BZYY  =
     >  - D2G*Y *.5D0 +4.D0*D4G*Y3/76.8D0  - 6.D0*D6G*Y5*Q7
      D3BX(3,1,1)= D3G*Y - D5G*Y3/12.D0
      D3BX(3,2,1)= D2G   - D4G*Y2*.25D0 +   D6G*Y4/76.8D0
      D3BX(3,2,2)=       - D3G*Y *.5D0  +4.D0*D5G*Y3/76.8D0
      D3BY(3,2,2)=  - D2G*.5D0 + D4G*Y2*.15625D0 - 30.D0*D6G*Y4*Q7

C     ... d2Bx/dXdY = d2By/dX2
      DDB(2,1,1)=
     > Z*D3BX(3,2,1) - D4G/12.D0*Z3 +  D6G/384.D0*(6.D0*Y2Z3 + Z5)
C     ... d2Bx/dY2
      DDB(2,2,1) = Z*D3BX(3,2,2) + D5G/32.D0*YZ3
C     ... d2Bx/dX2
      DDB(1,1,1) = Z*D3BX(3,1,1) - D5G/12.D0*YZ3
C     ... d2By/dY2
      DDB(2,2,2) =
     > Z*D3BY(3,2,2) + D4G/32.D0*Z3- D6G/128.D0*(Y2Z3 + 0.1D0*Z5)

C     .. dBx/dZ = dBz/dX
      DB(3,1)  = BZX - D3G/4.D0*YZ2 + D5G/76.8D0*(1.2D0*Y3Z2 + YZ4)
C     .. dBy/dZ = dBZ/dY
      DB(3,2)  = BZY - D2G/4.D0*Z2  + D4G/76.8D0*(3.6D0*Y2Z2 + Z4)
     >       + D6G/23040.D0*(45.D0*(Y2Z2+Z4) + 7.D0*Z6)

C     ... dBx/dX
      DB(1,1) = Z*BZXX - D4G/12.D0*YZ3 + D6G/384.D0*(2.D0*Y3Z3 + YZ5)
C     ... dBx/dY = dBy/dX
      DB(2,1) = Z*BZXY - D3G/12.D0*Z3  + D5G/384.D0*(6.D0*Y2Z3 + Z5)
C     ... dBy/dY
      DB(2,2) = Z*BZYY + D4G/32.D0*YZ3 - D6G/384.D0*(Y3Z3 + 0.3D0*YZ5)

C     .. d2By/dXdZ = d2Bz/dXdY = d2Bx/dYdZ
C   correction FM 7/12/94
C      DDB(3,2,1) = BZXY - D3G/4.D0*D0Z2  + D5G/76.8D0*(3.6D0*Y2Z2 + Z4)
      DDB(3,2,1) = BZXY - D3G/4.D0*Z2  + D5G/76.8D0*(3.6D0*Y2Z2 + Z4)
C     .. d2Bx/dXdZ = d2Bz/dX2
      DDB(3,1,1) = BZXX - D4G/4.D0*YZ2 + D6G/76.8D0*(1.2D0*Y3Z2 + YZ4)
C     .. d2By/dYdZ = d2Bz/dY2
      DDB(3,2,2) =
     > BZYY + D4G/32.D0*3.D0*YZ2 - D6G/256.D0*(2.D0*Y3Z2 + YZ4)

      IF(ID .GE. 3) THEN
        D4BX(3,1,1,1) =                D4G*Y       -    D6G*Y3/12.D0
        D4BX(3,2,1,1) = D3G             -   D5G*Y2*.25D0
        D4BX(3,2,2,1)=    -   D4G*Y*.5D0      + 4.D0*D6G*Y3/76.8D0
        D4BX(3,2,2,2)=    -   D3G  *.5D0     +12.D0*D5G*Y2/76.8D0
        D4BY(3,2,2,2)=         D4G*Y*.3125D0 - 14.D0*D6G*Y3/384.D0
        D4BX(3,3,3,1)=-(D4BX(3,1,1,1)  +D4BX(3,2,2,1))
        D4BX(3,3,3,2)=-(D4BX(3,2,1,1) +D4BX(3,2,2,2) )
        D4BY(3,3,3,2)=-(D4BX(3,2,2,1)+D4BY(3,2,2,2)  )
C       ... D3Bx/dX3
        D3BX(1,1,1) = Z*D4BX(3,1,1,1)   - D6G/12.D0 * YZ3
C       ... D3Bx/dX2DY
        D3BX(2,1,1) = Z*D4BX(3,2,1,1)  - D5G/12.D0 * Z3
C       ... D3Bx/dXdY2
        D3BX(2,2,1) = Z*D4BX(3,2,2,1) + D4G/32.D0 * YZ3
C       ... D3Bx/dY3
        D3BX(2,2,2) = Z*D4BX(3,2,2,2)  + D5G/32.D0 * Z3
C       ... D3By/dY3
        D3BY(2,2,2) = Z*D4BY(3,2,2,2)   - D6G/64.D0 * YZ3

      ENDIF
C     ... ENDIF ID>=3

      ENDIF

      ENDIF
C     **  TEST DLE+DLS

      IF(RT(IP) .NE. ZERO)
     >  CALL XROTB(RT(IP),B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
C-----  END QUAD

      CALL ADPOL(ID,IOP,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)

 83   CONTINUE
      IF( BM(3) .EQ. ZERO) GOTO 84

C--------------------------------------------------------------------
C----- Champ SEXTUPOLAIRE
      IP=3
      GN = BM(IP)*BRI
      DLE = DLEM(IP)
      DLS = DLSM(IP)
      Y=Y0
      Z=Z0
      IF(RT(IP) .NE. ZERO) CALL ROTX(RT(IP),Y,Z)

      IF(DLE+DLS .EQ. ZERO) THEN
C       *** SEXTUPOLE CRENEAU : PAS DE Champ DE FUITE

C        IF ( X.GE.XE .AND. X.LT.XS ) THEN
        IF ( X.GE.XE .AND. X.LE.XS ) THEN
C         *** A L'INTERIEUR DE L'AIMANT

          BZYY= 2.D0*GN
          DDB(3,2,2) = BZYY

          BZY = BZYY*Y
          DB(3,2) = BZY
          DB(2,2)= Z * BZYY

          B(1,2) = BZY*Z
          B(1,3) = GN*(Y+Z)*(Y-Z)

        ENDIF

      ELSEIF(DLE+DLS .NE. ZERO) THEN
C     *** SEXTUPOLE AVEC ChampS DE FUITE D'ENTREE ET/OU DE SORTIE

C     DISTANCE A LA FACE D'ENTREE :
      SE=XE-X
C     DISTANCE A LA FACE DE SORTIE :
      SS=X-XS

      ABSE=SE
      IF(ABSE.LT.ZERO) ABSE=-ABSE
      ABSS=SS
      IF(ABSS.LT.ZERO) ABSS=-ABSS

      CHFE=ABSE.LE.XE .AND. DLE .NE. ZERO
      CHFS=ABSS.LE.XLS .AND. DLS .NE. ZERO
      CHU=ABSE.GT.XE .AND. ABSS.GT.XLS

      IF(CHU) THEN
C       *** POINT DANS ZONE DE Champ Uniforme ( G(X)=0 )

          BZYY= 2.D0*GN
          DDB(3,2,2) = BZYY

          BZY = BZYY*Y
          DB(3,2) = BZY
          DB(2,2)= Z * BZYY

          B(1,2) = BZY*Z
          B(1,3) = GN*(Y+Z)*(Y-Z)
      ELSE
C     *** POINT DANS ZONE DE Champ DE FUITE D'ENTRE ET/OU DE SORTIE

      IF(.NOT. CHFE) THEN
C       *** HORS D'EFFET DU Champ DE FUITE D'ENTREE

        GE  =UN
        DGE =ZERO
        D2GE=ZERO
        D3GE=ZERO
        D4GE=ZERO

      ELSE
C       *** EFFET DU Champ DE FUITE D'ENTREE

        SE=SE/DLE
        CALL DRVG(4,CE,SE,GE,DGE,D2GE,D3GE,D4GE,D5GE,D6GE)

      ENDIF

      IF(.NOT. CHFS) THEN
C       *** HORS D'EFFET DU Champ DE FUITE DE SORTIE

        GS  =UN
        DGS =ZERO
        D2GS=ZERO
        D3GS=ZERO
        D4GS=ZERO

      ELSE
C       *** EFFET DU Champ DE FUITE DE SORTIE

        SS=SS/DLS
        CALL DRVG(4,CS,SS,GS,DGS,D2GS,D3GS,D4GS,D5GS,D6GS)

      ENDIF

      G= GE+GS-UN

      IF(G.LT.ZERO) CALL ENDJOB(
     >  'SBR CHAMC : gradient  is  wrong,  G (=GE+GS-1) < ',0)

C     *** CACULE LE GRADIENT DE Champ DE FUITE RESULTANT G=GE+GS-1,
C         ET  SES DERIVEES
      G   =  G*GN
      DG  = (DGE  * DE(IP,1) + DGS  * DS(IP,1))  *BRI
      D2G = (D2GE * DE(IP,2) + D2GS * DS(IP,2))  *BRI
      D3G = (D3GE * DE(IP,3) + D3GS * DS(IP,3))  *BRI
      D4G = (D4GE * DE(IP,4) + D4GS * DS(IP,4))  *BRI

      Y2 = Y  * Y
      Y3 = Y2 * Y
      Y4 = Y3 * Y

      Z2 = Z  * Z
      Z3 = Z2 * Z
      Z4 = Z3 * Z

      YZ   = Y  * Z
      YZ2  = Y  * Z2
      YZ3  = Y  * Z3

      Y2Z2 = Y2 * Z2

C----- FM, 03/93, ERROR CORRECTION
C      D4G=D4G*60.

C     *** Bx, By, Bz field components
      B(1,1) =Z*(DG*(Y2-Z2/3.D0) - D3G/48.D0*(3.D0*Y4 + 2.D0*Y2Z2 - Z4))
      B(1,2)  = YZ*(2.D0*G - D2G/12.D0*(3.D0*Y2 + Z2)
     >     + D4G/960.D0* (9.D0*Y4 + 10.D0*Y2Z2 + Z4))
      B(1,3)  =  G*(Y2 - Z2) - D2G/48.D0*(.6D0*Y4 + 1.2D0*Y2Z2 - Z4)
     >     + D4G/384.D0*(.6D0*Y4*Y2 + (3.D0*Y4 + Y2Z2 - 1.2D0*Z4)*Z2 )

C     *** DERIVEES DANS LE PLAN MEDIAN
      BZX   =   DG*Y2 - D3G*Y4 *.0625D0
      BZY   =  2.D0*G*Y - D2G*Y3 *.25D0
      BZXX  =  D2G*Y2 - D4G*Y4 *.0625D0
      BZXY  = 2.D0*DG*Y - D3G*Y3 *.25D0
      BZYY  =    2.D0*G - D2G*Y2 *.75D0
      D3BX(3,1,1) =  D3G*Y2
      D3BX(3,2,1) = 2.D0*D2G*Y - D4G*Y3 *.25D0
      D3BX(3,2,2) = 2.D0*DG    - D3G*Y2 *.75D0
      D3BY(3,2,2) =  - D2G*Y *1.5D0

C     ... d2Bx/dXdY = d2By/dX2
      DDB(2,1,1) = Z*D3BX(3,2,1) - D4G/12.D0*YZ3
C     ... d2Bx/dY2
      DDB(2,2,1) = Z*D3BX(3,2,2) - D3G/12.D0*Z3
C     ... d2Bx/dX2
      DDB(1,1,1) = Z*D3BX(3,1,1) - D3G/3.D0 *Z3
C     ... d2By/dY2
      DDB(2,2,2) = Z*D3BY(3,2,2)

C     .. dBx/dZ = dBz/dX
C      DB(3,1)  = BZX - DG*Z2 + D3G/9.6*(1.4*Y2Z2 - Z4)
      DB(3,1)  = BZX - DG*Z2 - D3G/9.6D0*(1.4D0*Y2Z3 - Z4)
C     .. dBy/dZ = dBz/dY
      DB(3,2)  = BZY - D2G*.25D0*YZ2

C     ... dBx/dX
C      DB(1,1) = Z*BZXX - D2G/3.D0*Z3 -D4G/48.D0*(2.D0*Y2Z2-Z*Z4)
C   correction FM 26/3/93
      DB(1,1) = Z*(BZXX - D2G/3.D0*Z2 -D4G/48.D0*(2.D0*Y2Z2-Z4))
C     ... dBx/dY = dBy/dX
      DB(2,1) = Z*BZXY - D3G/12.D0*YZ3
C     ... dBy/dY
      DB(2,2) = Z*(BZYY -D2G/12.D0*Z2+D4G/960.D0*(45.D0*Y4+
     > 30.D0*Y2Z2+Z4))

C     .. d2Bx/dXdZ = d2Bz/dX2
      DDB(3,1,1) = BZXX - D2G*Z2 - D4G*0.125D0*(Y2Z2 - Z4/1.2D0)
C     .. d2By/dXdZ = d2BZ/dXdY = d2Bx/dYdZ
      DDB(3,2,1) = BZXY - D3G*YZ2*.25D0
C     .. d2By/dYdZ = d2BZ/dY2
      DDB(3,2,2) = BZYY - .25D0*D2G*Z2

      IF(ID .GE. 3) THEN
        D4BX(3,1,1,1)  =    D4G*Y2
        D4BX(3,2,1,1) = 2.D0*D3G*Y
        D4BX(3,2,2,1)=   2.D0*D2G - D4G*Y2 *.75D0
        D4BX(3,2,2,2) =       - D3G*Y  *1.5D0
        D4BY(3,2,2,2)  =  -1.5D0*D2G
        D4BX(3,3,3,1)=-(D4BX(3,1,1,1)+  D4BX(3,2,2,1))
        D4BX(3,3,3,2)=-(D4BX(3,2,1,1)+ D4BX(3,2,2,2) )
        D4BY(3,3,3,2)=-(D4BX(3,2,2,1)+D4BY(3,2,2,2)  )

C       ... D3Bx/dX3
        D3BX(1,1,1) = Z*D4BX(3,1,1,1) - D4G * Z3/3D0
C       ... D3Bx/dX2DY
        D3BX(2,1,1) = Z*D4BX(3,2,1,1)
C       ... D3Bx/dXdY2
        D3BX(2,2,1) = Z*D4BX(3,2,2,1) - D4G *Z3/12.D0
C       ... D3Bx/dY3
        D3BX(2,2,2) = Z*D4BX(3,2,2,2)
C       ... D3By/dY3
        D3BY(2,2,2) = Z*D4BY(3,2,2,2)
      ENDIF
C     ... ENDIF ID>=3

      ENDIF

      ENDIF
C--------- TEST DLE+DLS

      IF(RT(IP) .NE. ZERO)
     >  CALL XROTB(RT(IP),B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
C------ END  SEXTU

      CALL ADPOL(ID,IOP,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)

 84   CONTINUE

      I3 = 3
      CALL ADPOL(ID,I3,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)

c                call zgnoel(
c     >            noel)
c                write(*,*) ' agsmmf bm ',noel,(bm(iu),iu=1,3)
c                write(*,*) ' agsmmf b(1,i) ',(b(1,iu),iu=1,3)
c                write(*,*) ' agsmmf gn, bri ',gn, bri, g
C                   read(*,*)
      RETURN
      END
