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
      SUBROUTINE MULTIP(ID,KFL,KUASEX,X,Y0,Z0,BM,DLEM,DLSM,DE,DS,RT
     >,XE,XS,CE,CS
     >,B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ,BT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(MPOL=10)
      DIMENSION BM(MPOL),DLEM(MPOL),DLSM(MPOL)
     >,DE(MPOL,6),DS(MPOL,6),RT(MPOL)
      PARAMETER(MCOEF=6)
      DIMENSION CE(MCOEF),CS(MCOEF)
      DIMENSION B(5,3),DB(3,3),DDB(3,3,3)
      DIMENSION D3BX(3,3,3), D3BY(3,3,3), D3BZ(3,3,3)
      DIMENSION D4BX(3,3,3,3) ,D4BY(3,3,3,3) ,D4BZ(3,3,3,3)
      DIMENSION BT(5, *)
 
      COMMON/AIM/ BO,RO,FG,GF,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXTRA.H"
      COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,BMAX
     > ,YCH,ZCH
      COMMON/CONST2/ ZERO, UN
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      COMMON/RIGID/ BORO,DPREF,DP,BR

      LOGICAL CHFE, CHFS, CHU
      PARAMETER (Q7=3.038194444D-4)

      LOGICAL CASPI, LTEMP
      SAVE CASPI

      XLS=XLIM-XS 
      IOP=1

      GOTO  (1,2,3,4,5,6,7,8,9,10,11), KUASEX
      CALL ENDJOB('*** Error, SBR MULTIP -> KUASEX value',-99)

 11   CONTINUE
      IF( BM(1) .EQ. ZERO) GOTO 82
 
C----- DIPOLE
 1    CONTINUE
 
      IP=1
      GN = BM(IP)/BR
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
          DG  = (DGE  * DE(IP,1) + DGS  * DS(IP,1))  /BR
          D2G = (D2GE * DE(IP,2) + D2GS * DS(IP,2))  /BR
          D3G = (D3GE * DE(IP,3) + D3GS * DS(IP,3))  /BR
          D4G = (D4GE * DE(IP,4) + D4GS * DS(IP,4))  /BR
          D5G = (D5GE * DE(IP,5) + D5GS * DS(IP,5))  /BR
          D6G = (D6GE * DE(IP,6) + D6GS * DS(IP,6))  /BR

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
c              write(*,*) 'sbr multip ',b(1,3),G ,D2G/8.D0, D4G/192.D0,
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
 
C 999  CONTINUE
 
      IF(RT(IP) .NE. ZERO)
     >  CALL XROTB(RT(IP),B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
      IF(KUASEX .EQ. IP) RETURN
C-----  END DIP

      CALL ADPOL(ID,IOP,KFL,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)

 82   CONTINUE
      IF( BM(2) .EQ. ZERO) GOTO 83
 
C--------------------------------------------------------------------
C----- QUADRUPOLE
 2    CONTINUE
 
      IP=2
      GN = BM(IP)/BR
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
      DG  = (DGE  * DE(IP,1) + DGS  * DS(IP,1))  /BR
      D2G = (D2GE * DE(IP,2) + D2GS * DS(IP,2))  /BR
      D3G = (D3GE * DE(IP,3) + D3GS * DS(IP,3))  /BR
      D4G = (D4GE * DE(IP,4) + D4GS * DS(IP,4))  /BR
      D5G = (D5GE * DE(IP,5) + D5GS * DS(IP,5))  /BR
      D6G = (D6GE * DE(IP,6) + D6GS * DS(IP,6))  /BR
 
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
      IF(KUASEX .EQ. IP) RETURN
C-----  END QUAD
 
      CALL ADPOL(ID,IOP,KFL,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)

 83   CONTINUE
      IF( BM(3) .EQ. ZERO) GOTO 84
 
C--------------------------------------------------------------------
C----- Champ SEXTUPOLAIRE
 3    CONTINUE
 
      IP=3
      GN = BM(IP)/BR
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
      DG  = (DGE  * DE(IP,1) + DGS  * DS(IP,1))  /BR
      D2G = (D2GE * DE(IP,2) + D2GS * DS(IP,2))  /BR
      D3G = (D3GE * DE(IP,3) + D3GS * DS(IP,3))  /BR
      D4G = (D4GE * DE(IP,4) + D4GS * DS(IP,4))  /BR
 
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
      IF(KUASEX .EQ. IP) RETURN
C------ END  SEXTU
 
      CALL ADPOL(ID,IOP,KFL,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)
 
 84   CONTINUE
      IF( BM(4) .EQ. ZERO) GOTO 85
 
C--------------------------------------------------------------------
C----- OCTUPOLE
 4    CONTINUE
 
      IP=4
      GN = BM(IP)/BR
      DLE = DLEM(IP)
      DLS = DLSM(IP)
      Y=Y0
      Z=Z0
      IF(RT(IP) .NE. ZERO) CALL ROTX(RT(IP),Y,Z)
 
      IF(DLE+DLS .EQ. ZERO) THEN
C       *** OCTUPOLE CRENEAU : PAS DE Champ DE FUITE
 
C        IF ( X.GE.XE .AND. X.LT.XS ) THEN
        IF ( X.GE.XE .AND. X.LE.XS ) THEN
C         *** A L'INTERIEUR DE L'AIMANT
 
          D3BY(3,2,2) = 6.D0*GN
 
          BZYY= Y * D3BY(3,2,2)
          DDB(3,2,2) = BZYY
          DDB(2,2,2) = Z * D3BY(3,2,2)
 
          BZY = BZYY*Y*.5D0
          DB(2,2) = Z * BZYY
          DB(3,2) = BZY - 3.D0*Z*Z*GN
 
          B(1,2) = Z *(BZY - GN*Z*Z)
          B(1,3) = GN * Y * (Y*Y - 3.D0*Z*Z)
        ENDIF
 
      ELSEIF(DLE+DLS .NE. ZERO) THEN
C     *** OCTUPOLE AVEC ChampS DE FUITE D'ENTREE ET/OU DE SORTIE
 
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
      CHU= ABSE.GT.XE .AND. ABSS.GT.XLS
 
      IF(CHU) THEN
C       *** POINT DANS ZONE DE Champ Uniforme ( G(X)=0 )
 
          D3BY(3,2,2) = 6.D0*GN
C          D3BY(3,3,3)  = - BYYY
 
          BZYY= Y * D3BY(3,2,2)
          DDB(3,2,2) = BZYY
          DDB(2,2,2) = Z * D3BY(3,2,2)
 
          BZY = BZYY*Y*.5D0
          DB(2,2) = Z * BZYY
          DB(3,2) = BZY - 3.D0*Z*Z*GN
 
          B(1,2) = Z *(BZY - GN*Z*Z)
          B(1,3) = GN * Y * (Y*Y - 3.D0*Z*Z)
      ELSE
C     *** POINT DANS ZONE DE Champ DE FUITE D'ENTRE ET/OU DE SORTIE
 
      IF(.NOT. CHFE) THEN
C       *** HORS D'EFFET DU Champ DE FUITE D'ENTREE
C           OU GRADIENT CRENEAU
 
        GE  =UN
        DGE =ZERO
        D2GE=ZERO
        D3GE=ZERO
        D4GE=ZERO
 
      ELSE
C       *** EFFET DU Champ DE FUITE D'ENTREE
 
        SE=SE/DLE
        CALL DRVG(6,CE,SE,GE,DGE,D2GE,D3GE,D4GE,D5GE,D6GE)
 
      ENDIF
 
      IF(.NOT. CHFS) THEN
C       *** HORS D'EFFET DU Champ DE FUITE DE SORTIE
C           OU GRADIENT CRENEAU
 
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
     >  'SBR CHAMC : gradient is  wrong,  G (=GE+GS-1) < ',0)
 
C     *** CACULE LE DE Champ DE FUITE RESULTANT G=GE+GS-1, ET
C         SES DERIVEES
      G   =  G*GN
      DG  = (DGE  * DE(IP,1) + DGS  * DS(IP,1))  /BR
      D2G = (D2GE * DE(IP,2) + D2GS * DS(IP,2))  /BR
      D3G = (D3GE * DE(IP,3) + D3GS * DS(IP,3))  /BR
      D4G = (D4GE * DE(IP,4) + D4GS * DS(IP,4))  /BR
 
      Y2 = Y  * Y
      Y3 = Y2 * Y
      Y4 = Y3 * Y
      Y5 = Y4 * Y
 
      Z2 = Z  * Z
      Z3 = Z2 * Z
      Z4 = Z3 * Z
 
      YZ   = Y  * Z
      YZ2  = Y  * Z2
      YZ3  = Y  * Z3
      YZ4  = Y  * Z4
 
      Y2Z2 = Y2 * Z2
 
C     *** COMPOSANTES Bx, By, Bz DU Champ
      B(1,1)  =  YZ*(DG*(Y2-Z2) - D3G*.05D0*(Y4 - Z4))
      B(1,2)  =  Z*(G*(3.D0*Y2 - Z2) - D2G*.05D0*(5.D0*Y4 - Z4)
     >     + D4G/960.D0*( (7.D0*Y4 + 5.D0*Y2Z2 - 3.D0*Z4)*Y2 - Z3*Z3 ))
      B(1,3)  =  G*(Y3 - 3.D0*YZ2) - D2G*.05D0*(Y5 - YZ4)
     >     + D4G/960.D0*( Y4*Y3 + (3.D0*Y4 - 5.D0*Y2Z2 - 7.D0*Z4) *YZ2 )
 
C     *** DERIVEES DANS LE PLAN MEDIAN
      BZX   =   DG*Y3 - D3G*Y5 *.05D0
      BZY   =  3.D0*G*Y2- D2G*Y4 *.25D0
      BZXX  =  D2G*Y3 - D4G*Y5 *.05D0
      BZXY  = 3.D0*DG*Y2- D3G*Y4 *.25D0
      BZYY  = 6.D0* G*Y - D2G*Y3
      D3BX(3,1,1) =  D3G*Y3
      D3BX(3,2,1) = 3.D0*D2G*Y2- D4G*Y4 *.25D0
      D3BX(3,2,2) = 6.D0*DG *Y - D3G*Y3
      D3BY(3,2,2) = 6.D0* G    - D2G*Y2 *3D0
 
C     ... d2Bx/dXdY = d2By/dX2
      DDB(2,1,1) = Z*D3BX(3,2,1) - D2G*Z3
C     ... d2Bx/dY2
      DDB(2,2,1) = Z*D3BX(3,2,2)
C     ... d2Bx/dX2
      DDB(1,1,1) = Z*D3BX(3,1,1) - D3G *YZ3
C     ... d2By/dY2
      DDB(2,2,2) = Z*D3BY(3,2,2)
 
C     .. dBx/dZ = dBz/dX
      DB(3,1)  = BZX - 3.D0*DG*YZ2 + D3G*.25D0*YZ4
C     .. dBy/dZ = dBz/dY
      DB(3,2)  = BZY - 3.D0*G*Z2   + 1.25D0*D2G*Z4
 
C     ... dBx/dX
      DB(1,1) = Z*BZXX - D2G*YZ3
C     ... dBx/dY = dBy/dX
      DB(2,1) = Z*(BZXY  + D3G*Z4*.05D0) - DG*Z3
C     ... dBy/dY
      DB(2,2) = Z*BZYY
 
C     .. d2Bx/dXdZ = d2Bz/dX2
      DDB(3,1,1) = BZXX - 3.D0*D2G*YZ2 + .25D0*D4G*YZ4
C     .. d2By/dXdZ = d2Bz/dXdY = d2Bx/dYdZ
      DDB(3,2,1) = BZXY - 3.D0*DG*Z2 + D3G*Z4*.25D0
C     .. d2By/dYdZ = d2Bz/dY2
      DDB(3,2,2) = BZYY

      IF(ID .GE. 3) THEN
        D4BX(3,1,1,1)  =    D4G*Y3
        D4BX(3,2,1,1) = 3.D0*D3G*Y2
        D4BX(3,2,2,1)=   6.D0*D2G*Y - D4G*Y3
        D4BX(3,2,2,2) = 6.D0*DG - 3.D0*D3G*Y2
        D4BY(3,2,2,2)  =  -6.D0*D2G*Y
        D4BX(3,3,3,1)=-(D4BX(3,1,1,1)+D4BX(3,2,2,1))
        D4BX(3,3,3,2)=-(D4BX(3,2,1,1)+D4BX(3,2,2,2))
        D4BY(3,3,3,2)=-(D4BX(3,2,2,1)+D4BY(3,2,2,2))
 
C       ... D3Bx/dX3
        D3BX(1,1,1) = Z*D4BX(3,1,1,1) - D4G*YZ3
C       ... D3Bx/dX2DY
        D3BX(2,1,1) = Z*D4BX(3,2,1,1)
C       ... D3Bx/dXdY2
        D3BX(2,2,1) = Z*D4BX(3,2,2,1)
C       ... D3Bx/dY3
        D3BX(2,2,2) = Z*D4BX(3,2,2,2)
C       ... D3By/dY3
        D3BY(2,2,2) = Z*D4BY(3,2,2,2)
      ENDIF
C     ... ENDIF ID>=3
 
      ENDIF
 
      ENDIF
C     ** TEST DLE+DLS
 
      IF(RT(IP) .NE. ZERO)
     >  CALL XROTB(RT(IP),B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
      IF(KUASEX .EQ. IP) RETURN
C-----  END OCTUPOLE
 
      CALL ADPOL(ID,IOP,KFL,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)
 
 85   CONTINUE
      IF( BM(5) .EQ. ZERO) GOTO 86
 
C--------------------------------------------------------------------
C----- DECAPOLE
 5    CONTINUE
 
      IP=5
      GN = BM(IP)/BR
      DLE = DLEM(IP)
      DLS = DLSM(IP)
      Y=Y0
      Z=Z0
      IF(RT(IP) .NE. ZERO) CALL ROTX(RT(IP),Y,Z)
 
      IF(DLE+DLS .EQ. ZERO) THEN
C       *** DECAPOLE CRENEAU : PAS DE Champ DE FUITE
 
C        IF ( X.GE.XE .AND. X.LT.XS ) THEN
        IF ( X.GE.XE .AND. X.LE.XS ) THEN
C         *** A L'INTERIEUR DE L'AIMANT
 
          D4BY(3,2,2,2) = 24.D0*GN
          D3BY(3,2,2) = Y*D4BY(3,2,2,2)
          BZYY = .5D0*Y*D3BY(3,2,2)
          BZY = Y*BZYY/3D0
 
          DDB(2,2,2) = Z * D3BY(3,2,2)
 
          Z2=Z*Z
          DB(3,2) = BZY - 12.D0*Y*Z2*GN
          DB(2,2) = Z * (BZYY - 4.D0*GN*Z2)
 
          DDB(3,2,2) = BZYY - 12.D0*GN*Z2
 
          B(1,2) = Z *(BZY - 4.D0*GN*Y*Z2)
          B(1,3) = GN*(Y*Y*(Y*Y-6.D0*Z2) + Z2*Z2)

          IF(ID.GE.3) THEN
            D4BY(3,3,3,2)=-D4BY(3,2,2,2)
            D3BY(2,2,2) = Z*D4BY(3,2,2,2)
          ENDIF

        ENDIF
 
      ELSEIF(DLE+DLS .NE. ZERO) THEN
C     *** DECAPOLE AVEC ChampS DE FUITE D'ENTREE ET/OU DE SORTIE
 
        CALL ENDJOB('Decapole fringe field not installed ',-99)
 
      ENDIF
C     ** TEST DLE+DLS
 
      IF(RT(IP) .NE. ZERO)
     >  CALL XROTB(RT(IP),B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
      IF(KUASEX .EQ. IP) RETURN
C-----  END DECAPOLE
 
      CALL ADPOL(ID,IOP,KFL,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)
 
 86   CONTINUE
      IF( BM(6) .EQ. ZERO) GOTO 87
 
C--------------------------------------------------------------------
C----- DODECAPOLE
 6    CONTINUE
 
      IP=6
      GN = BM(IP)/BR
      DLE = DLEM(IP)
      DLS = DLSM(IP)
      Y=Y0
      Z=Z0
      IF(RT(IP) .NE. ZERO) CALL ROTX(RT(IP),Y,Z)
 
      IF(DLE+DLS .EQ. ZERO) THEN
C       *** DODECAPOLE CRENEAU 
 
C        IF ( X.GE.XE .AND. X.LT.XS ) THEN
        IF ( X.GE.XE .AND. X.LE.XS ) THEN
C         *** A L'INTERIEUR DE L'AIMANT
 
          D4BY(3,2,2,2) = 120.D0*GN*Y
          D3BY(3,2,2) = .5D0*Y*D4BY(3,2,2,2)
          BZYY = Y*D3BY(3,2,2)/3D0
          BZY = .25D0*Y*BZYY
 
          Z2=Z*Z
          SGZ2=60.D0*GN*Z2
 
          DDB(2,2,2) = Z * D3BY(3,2,2) - 20.D0*GN*Z2*Z
 
          Y2=Y*Y
          DB(3,2) = BZY + 5.D0*Z2*GN*(-6.D0*Y2+Z2)
          DB(2,2) = 20.D0*GN*Y*Z*(Y2-Z2)
 
          DDB(3,2,2) = BZYY - 60.D0*GN*Y*Z2
 
          B(1,2) = Z *(BZY + (-10.D0*Y2 + Z2)*GN*Z2)
          B(1,3) = GN*Y*(Y2*(Y2-10.D0*Z2) + 5.D0*Z2*Z2 )

          IF(ID.GE.3) THEN
            D4BY(3,3,3,2)=-D4BY(3,2,2,2)
            D3BY(2,2,2) = Z*D4BY(3,2,2,2)
          ENDIF

        ENDIF
 
      ELSEIF(DLE+DLS .NE. ZERO) THEN
C     *** DODECAPOLE AVEC ChampS DE FUITE D'ENTREE ET/OU DE SORTIE
 
        CALL ENDJOB('Dodecapole fringe field not installed ',-99)
 
      ENDIF
C     ** TEST DLE+DLS
 
      IF(RT(IP) .NE. ZERO)
     >  CALL XROTB(RT(IP),B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
      IF(KUASEX .EQ. IP) RETURN
C-----  END DODECAPOLE
 
      CALL ADPOL(ID,IOP,KFL,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)

 87   CONTINUE
      IF( BM(7) .EQ. ZERO) GOTO 88
 
C--------------------------------------------------------------------
C----- 14-POLE
 7    CONTINUE
 
      IP=7
      GN = BM(IP)/BR
      DLE = DLEM(IP)
      DLS = DLSM(IP)
      Y=Y0
      Z=Z0
      IF(RT(IP) .NE. ZERO) CALL ROTX(RT(IP),Y,Z)
 
      IF(DLE+DLS .EQ. ZERO) THEN
C       *** 14-POLE CRENEAU 
 
C        IF ( X.GE.XE .AND. X.LT.XS ) THEN
        IF ( X.GE.XE .AND. X.LE.XS ) THEN
C         *** A L'INTERIEUR DE L'AIMANT

          CALL ENDJOB('14-pole not installed ',-99)
 
          IF(ID.GE.3) THEN
          ENDIF

        ENDIF
 
      ELSEIF(DLE+DLS .NE. ZERO) THEN
C     *** 14-POLE AVEC ChampS DE FUITE D'ENTREE ET/OU DE SORTIE
 
        CALL ENDJOB('14-pole fringe field not installed ',-99)
 
      ENDIF
C     ** TEST DLE+DLS
 
      IF(RT(IP) .NE. ZERO)
     >  CALL XROTB(RT(IP),B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
      IF(KUASEX .EQ. IP) RETURN
C-----  END 14-POLE
 
      CALL ADPOL(ID,IOP,KFL,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)

 88   CONTINUE
      IF( BM(8) .EQ. ZERO) GOTO 89
 
C--------------------------------------------------------------------
C----- 16-POLE
 8    CONTINUE
 
      IP=8
      GN = BM(IP)/BR
      IF(GN.NE.0.D0) 
     >  stop ' *** 16-POLE lens not installed, please set B7=0 ***     '

      DLE = DLEM(IP)
      DLS = DLSM(IP)
      Y=Y0
      Z=Z0
      IF(RT(IP) .NE. ZERO) CALL ROTX(RT(IP),Y,Z)
 
      IF(DLE+DLS .EQ. ZERO) THEN
C       *** CRENEAU : PAS DE Champ DE FUITE
 
C        IF ( X.GE.XE .AND. X.LT.XS ) THEN
        IF ( X.GE.XE .AND. X.LE.XS ) THEN
C         *** A L'INTERIEUR DE L'AIMANT
 
        ENDIF
 
      ELSEIF(DLE+DLS .NE. ZERO) THEN
C     *** ChampS DE FUITE D'ENTREE ET/OU DE SORTIE
 
        CALL ENDJOB('SBR MULTIP:  fringe field not installed ',-99)
 
      ENDIF
C     ** TEST DLE+DLS
 
      IF(RT(IP) .NE. ZERO)
     >  CALL XROTB(RT(IP),B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
      IF(KUASEX .EQ. 5) RETURN
C-----  END 16-POLE
 
      CALL ADPOL(ID,IOP,KFL,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)
 
 89   CONTINUE
      IF( BM(9) .EQ. ZERO) GOTO 90
 
C--------------------------------------------------------------------
C----- 18-POLE
 9    CONTINUE
 
      IP=9
      GN = BM(IP)/BR
      IF(GN.NE.0.D0) 
     >  stop ' *** 18-POLE lens not installed, please set B7=0 ***'

      DLE = DLEM(IP)
      DLS = DLSM(IP)
      Y=Y0
      Z=Z0
      IF(RT(IP) .NE. ZERO) CALL ROTX(RT(IP),Y,Z)
 
      IF(DLE+DLS .EQ. ZERO) THEN
C       *** CRENEAU : PAS DE Champ DE FUITE
 
C        IF ( X.GE.XE .AND. X.LT.XS ) THEN
        IF ( X.GE.XE .AND. X.LE.XS ) THEN
C         *** A L'INTERIEUR DE L'AIMANT
 
        ENDIF
 
      ELSEIF(DLE+DLS .NE. ZERO) THEN
C     *** ChampS DE FUITE D'ENTREE ET/OU DE SORTIE
 
        CALL ENDJOB('SBR MULTIP:  fringe field not installed ',-99)
 
      ENDIF
C     ** TEST DLE+DLS
 
      IF(RT(IP) .NE. ZERO)
     >  CALL XROTB(RT(IP),B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
      IF(KUASEX .EQ. 5) RETURN
C-----  END DECAPOLE
 
      CALL ADPOL(ID,IOP,KFL,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)
 
 90   CONTINUE
      IF( BM(10) .EQ. ZERO) GOTO 91
 
C--------------------------------------------------------------------
C----- 20-POLE
 10   CONTINUE
 
      IP=10
      GN = BM(IP)/BR
      DLE = DLEM(IP)
      DLS = DLSM(IP)
      Y=Y0
      Z=Z0
      IF(RT(IP) .NE. ZERO) CALL ROTX(RT(IP),Y,Z)
 
      IF(DLE+DLS .EQ. ZERO) THEN
C       *** CRENEAU : PAS DE Champ DE FUITE
 
C        IF ( X.GE.XE .AND. X.LT.XS ) THEN
        IF ( X.GE.XE .AND. X.LE.XS ) THEN
C         *** A L'INTERIEUR DE L'AIMANT
 
          
          Y2=Y*Y
          Y4 = Y2*Y2
          Y6 = Y4*Y2
          Y8 = Y4*Y4
          Z2=Z*Z
          Z4 = Z2*Z2
          Z6 = Z4*Z2
          Z8 = Z4*Z4
          DDB(2,2,2)=GN*Z*(504.D0*Y6-420.D0*Y4*Z2+1512.D0*Y2*Z4-72D0*Z6)
 
          DB(3,2)=GN*(9.D0*Y8-252.D0*Y6*Z2+630.D0*Y4*Z4-252.D0*Y2*Z6+Z8)
          DB(2,2)=GN*Z*Y*(72.D0*Y6-84.D0*Y4*Z2+504.D0*Y2*Z4-72.D0*Z6)
 
          DDB(3,2,2)=GN*Y*(72.D0*Y6-1512.D0*Y4*Z2+2520D0*Y2*Z4-504D0*Z6)
 
          B(1,2) =GN*Z*(9.D0*Y8-84.D0*Y6*Z2+126.D0*Y4*Z4-36.D0*Y2*Z6+Z8)
          B(1,3) =GN*Y*(Y8-36.D0*Y6*Z2+126.D0*Y4*Z4-84.D0*Y2*Z6+9.D0*Z8)

        ENDIF
 
      ELSEIF(DLE+DLS .NE. ZERO) THEN
C     *** ChampS DE FUITE D'ENTREE ET/OU DE SORTIE
 
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
      CHU= ABSE.GT.XE .AND. ABSS.GT.XLS
 
      IF(CHU) THEN
C       *** POINT DANS ZONE DE Champ Uniforme ( G(X)=0 )
 
        IF(.NOT. CASPI) THEN

          Y2=Y*Y
          Y4 = Y2*Y2
          Y6 = Y4*Y2
          Y8 = Y4*Y4
          Z2=Z*Z
          Z4 = Z2*Z2
          Z6 = Z4*Z2
          Z8 = Z4*Z4
          DDB(2,2,2)=GN*Z*(504.D0*Y6-420.D0*Y4*Z2+1512.D0*Y2*Z4-72D0*Z6)
 
          DB(3,2)=GN*(9.D0*Y8-252.D0*Y6*Z2+630.D0*Y4*Z4-252.D0*Y2*Z6+Z8)
          DB(2,2)=GN*Z*Y*(72.D0*Y6-84.D0*Y4*Z2+504.D0*Y2*Z4-72.D0*Z6)
 
          DDB(3,2,2)=GN*Y*(72.D0*Y6-1512.D0*Y4*Z2+2520D0*Y2*Z4-504D0*Z6)
 
          B(1,2) =GN*Z*(9.D0*Y8-84.D0*Y6*Z2+126.D0*Y4*Z4-36.D0*Y2*Z6+Z8)
          B(1,3) =GN*Y*(Y8-36.D0*Y6*Z2+126.D0*Y4*Z4-84.D0*Y2*Z6+9.D0*Z8)

        ELSE
          DDB(2,2,2)=0
 
          DB(3,2)=GN*0.D0
          DB(2,2)=GN*0.D0
          DDB(3,2,2)=0.D0
          B(1,2) =GN*0.D0
          B(1,3) =GN*0.D0

        ENDIF

      ELSE
C     *** POINT DANS ZONE DE Champ DE FUITE D'ENTRE ET/OU DE SORTIE
 


      IF(.NOT. CHFE) THEN
C       *** HORS D'EFFET DU Champ DE FUITE D'ENTREE
C           OU GRADIENT CRENEAU
 
        GE  =UN
        DGE =ZERO
        D2GE=ZERO
        D3GE=ZERO
 
      ELSE
C       *** EFFET DU Champ DE FUITE D'ENTREE
 
        SE=SE/DLE
        CALL DRVG(6,CE,SE,GE,DGE,D2GE,D3GE,D4GE,D5GE,D6GE)
 
      ENDIF
 
      IF(.NOT. CHFS) THEN
C       *** HORS D'EFFET DU Champ DE FUITE DE SORTIE
C           OU GRADIENT CRENEAU
 
        GS  =UN
        DGS =ZERO
        D2GS=ZERO
        D3GS=ZERO
 
      ELSE
C       *** EFFET DU Champ DE FUITE DE SORTIE
 
        SS=SS/DLS
        CALL DRVG(4,CS,SS,GS,DGS,D2GS,D3GS,D4GS,D5GS,D6GS)
 
      ENDIF
 
      G= GE+GS-UN
 
      IF(G.LT.ZERO) CALL ENDJOB(
     >  'SBR CHAMC : GRADIENT ERRONNE  G (=GE+GS-1) < ',0)
 
C     *** CACULE LE DE Champ DE FUITE RESULTANT G=GE+GS-1, ET
C         SES DERIVEES
      G   =  G*GN
      DG  = (DGE  * DE(IP,1) + DGS  * DS(IP,1))  /BR
      D2G = (D2GE * DE(IP,2) + D2GS * DS(IP,2))  /BR
      D3G = (D3GE * DE(IP,3) + D3GS * DS(IP,3))  /BR
 
      G   =  G  /10.D0
      DG  = DG  /10.D0
      D2G = D2G /10.D0
      D3G = D3G /10.D0
 
      Y2 = Y  * Y
      Y4 = Y2 * Y2
      Y6 = Y4 * Y2
      Y8 = Y4 * Y4
 
      Z2 = Z  * Z
      Z4 = Z2 * Z2
      Z6 = Z4 * Z2
      Z8 = Z4 * Z4
 
      YZ   = Y  * Z
 
      IF(.NOT. CASPI) THEN

      U = G - D2G/44.D0
      V = 10.D0*Y8-120.D0*Y6*Z2+252.D0*Y4*Z4-120.D0*Y2*Z6+10.D0*Z8
      DUX = G - D3G/44.D0*(Y2+Z2)
      DUY = -D2G/22.D0*Y
      DUZ = -D2G/22.D0*Z
C     *** COMPOSANTES Bx, By, Bz DU Champ
      B(1,1)  =  DUX*V*YZ
      B(1,2)  =  DUY*V*YZ + U*DVY*YZ + U*V*Z
      B(1,3)  =  DUZ*V*YZ + U*DVZ*YZ + U*V*Y
 
C     Derivees
      DUXX = D2G
      DUXY = -D3G/22.D0*Y
      DUXZ = -D3G/22.D0*Z
      DUYY = -D2G/22.D0
      DUZZ = -D2G/22.D0
      DUYZ = 0.D0
      DUXXX = D3G
      DUXYY = -D3G/22.D0
      DVY = (80.D0*Y6-720.D0*Y4*Z2+1008.D0*Y2*Z4-240.D0*Z6)*Y
      DVZ = (80.D0*Z6-720.D0*Z4*Y2+1008.D0*Z2*Y4-240.D0*Y6)*Z
      DVYY = 560.D0*Y6 -3600.D0*Y4*Z2 + 3024.D0*Y2*Z4 - 240.D0*Z6
      DVZZ = 560.D0*Z6 -3600.D0*Z4*Y2 + 3024.D0*Z2*Y4 - 240.D0*Y6
      DVYZ = (-1440.D0*Y4 + 4032.D0*Y2*Z2 - 1440.D0*Z4)*YZ
      DVYYZ = (-7200.D0*Y4 + 12096.D0*Y2*Z2 - 1440.D0*Z4)*Z
      DVYYY = (3360.D0*Y4 - 14400.D0*Y2*Z2 + 6048.D0*Z4)*Y

C     ... d2Bx/dXdY = d2By/dX2
      DDB(2,1,1) = DUXX*(DVY*YZ + V*Z)
C     ... d2Bx/dY2
      DDB(2,2,1) = DUXYY*V*YZ + 2D0*DUXY*DVY*YZ + 2D0*DUXY*V*Z +
     >   DUX*DVYY*YZ
C     ... d2Bx/dX2
      DDB(1,1,1) = DUXXX*V*YZ
C     ... d2By/dY2
      DDB(2,2,2) = DUYY*V*Z + 3D0*DUYY*DVY*YZ + 6D0*DUY*DVY*Z +
     >   2D0*DUYY*V*Z + 3D0*DUY*DVYY*YZ + U*DVYYY*YZ + 3D0*U*DVYY*Z
 
C     .. dBx/dZ = dBz/dX
      DB(3,1)  = DUXZ*V*YZ + DUX*DVZ*YZ + DUX*V*Y
C     .. dBy/dZ = dBz/dY
      DB(3,2)  = DUY*DVZ*YZ + DUY*V*Y + DUZ*DVY*YZ+ U*DVYZ*YZ+ U*DVY*Y+
     >   DUZ*V*Z + U*DVZ*Z + U*V
C     ... dBx/dX
      DB(1,1) = DUXX*V*YZ
C     ... dBx/dY = dBy/dX
      DB(2,1) = DUXY*V*YZ + DUX*DVY*YZ + DUX*V*Z
C     ... dBy/dY
      DB(2,2) = DUYY*V*YZ + 2D0*DUY*DVY*YZ + 2D0*DUY*V*Z + U*DVYY*YZ +
     >  2D0*U*DVY*Z
 
C     .. d2Bx/dXdZ = d2Bz/dX2
      DDB(3,1,1) = DUXX*(DVZ*YZ + V*Y)
C     .. d2By/dXdZ = d2Bz/dXdY = d2Bx/dYdZ
      DDB(3,2,1) = DUXZ*DVY*YZ + DUXZ*V*Z + DUXY*DVZ*YZ + DUX*DVYZ*YZ +
     >   DUX*DVZ*Z + DUXY*V*Y + DUX*DVY*Y + DUX*V
C     .. d2By/dYdZ = d2Bz/dY2
      DDB(3,2,2) = DUYY*DVZ*YZ + DUYY*V*Y + 2D0*DUY*(DVYZ*YZ + DVY*Y + 
     >  DVZ*Z + V) + DUZ*DVYY*YZ + U*(DVYYZ*YZ + DVYY*Y + 2D0*(DVYZ*Z
     >  +DVY)) + 2D0*DUZ*DVY*Z
 

      ELSEIF(CASPI) THEN
C------------ For DA studies on LHC low-beta quads
        GN = GN / (2.D0*(XE + XLS)) * (XS-XE)
        B(1,2) =GN*Z*(9.D0*Y8-84.D0*Y6*Z2+126.D0*Y4*Z4-36.D0*Y2*Z6+Z8)
        B(1,3) =GN*Y*(Y8-36.D0*Y6*Z2+126.D0*Y4*Z4-84.D0*Y2*Z6+9.D0*Z8)

      ENDIF
C---------- CASPI


      ENDIF
C---------- CHU


      ENDIF
C---------- TEST DLE+DLS
 
      IF(RT(IP) .NE. ZERO)
     >  CALL XROTB(RT(IP),B,DB,DDB,D3BX,D3BY,D3BZ,D4BX,D4BY,D4BZ)
      IF(KUASEX .EQ. 5) RETURN
C-----  END 20-POLE
 
      CALL ADPOL(ID,IOP,KFL,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)
 
 91   CONTINUE
      I3 = 3
      CALL ADPOL(ID,I3,KFL,B,DB,DDB,D3BX,D3BY,D4BX,D4BY,BT)
 
      RETURN

      ENTRY MULTI1(LTEMP)
      CASPI = LTEMP
C      WRITE(6,*) ' CASPI SWITCHED ON !'
      RETURN

      END
