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
      SUBROUTINE ELCYL(IDE,E0,MPOL,QLEM,QLSM,QE,QS,RTQ,A,R, 
     >                                                     E,DE,DDE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QLEM(MPOL),QLSM(MPOL),QE(MPOL,6),QS(MPOL,6),RTQ(MPOL)
      DIMENSION E(5,3),DE(3,3),DDE(3,3,3)

      COMMON/AIM/ AE,AT,AS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER(MCOEF=6)
      COMMON/CHAFUI/ XE,XS,CE(MCOEF),CS(MCOEF),QCE(MCOEF),QCS(MCOEF)
      COMMON/INTEG/ PAS,DXI,XLIM,XCE,YCE,ALE,XCS,YCS,ALS,KP
      COMMON/RIGID/ BORO,DPREF,DP,BR

C------- Extent of Entrance/Exit fringe field region
      EQUIVALENCE (EB1,XLE), (EB2,XLS), (EG1,V0)

      LOGICAL CHFE, CHU, CHFS
      PARAMETER (I0=0)
 
      VN = V0/BR

C----- LambdaE,S
      QLE = QLEM(1)
      QLS = QLSM(1)

      IF(QLE+QLS .EQ. 0D0) THEN
C------- Sharp edges  !!! test !!!  non physical -> dp/p altered
 
        IF ( A .GE. AE .AND. A .LE. AE+AT ) THEN
C         Inside the deflector

          E(1,2) = - VN / R
C          DE(2,2) = - E(1,2) / R
C          DDE(2,2,2) = - DE(2,2) * 2. / R

        ENDIF

      ELSE
 
C        Position w.r.t. Entrance/Exit EFB 

        CHFE= A .LE. 2D0*AE .AND. QLE.NE.0D0
        CHU=  A .GT. 2D0*AE .AND. A .LT. AT - 2D0*AS   
        CHFS= A .GE. AT - 2D0*AS .AND. QLS.NE.0D0
  
        IF(CHU) THEN
C---------- Central region
 
          E(1,2) = - VN / R
C          DE(2,2) = - E(1,2) / R
C          DDE(2,2,2) = -2./R * DE(2,2) 
      
        ELSE
C---------- ENtrance or exit fringe field
 
          IF(.NOT. CHFE) THEN
C           *** HORS D'EFFET DU Champ DE FUITE D'ENTREE
C               OU GRADIENT CRENEAU
 
            GE  =1D0
            DGE = 0D0
            D2GE = 0D0

          ELSE
C           *** EFFET DU Champ DE FUITE D'ENTREE
 
C------------ SE>0 outside
            SE=(AE-A)*R/QLE
            CALL DRVG(IDE,QCE,SE,GE,DGE,D2GE,D3GE,D4GE,D5GE,D6GE)
C     >                                   D7GE,D8GE,D9GE,D10GE)

            DSEA = - 1.D0

          ENDIF
 
          IF(.NOT. CHFS) THEN
C           *** HORS D'EFFET DU Champ DE FUITE DE SORTIE
C               OU GRADIENT CRENEAU
  
            GS  =1D0
            DGS = 0D0
            D2GS = 0D0

          ELSE
C           *** EFFET DU Champ DE FUITE DE SORTIE
 
C------------ SS>0 outside
            SS=(A-(AT-AS))*R/QLS        
            CALL DRVG(IDE,QCS,SS,GS,DGS,D2GS,D3GS,D4GS,D5GS,D6GS)
C     >                                   D7GS,D8GS,D9GS,D10GS)

            DSSA = 1.D0

          ENDIF
 
          G = ( GE+GS-1D0 )

          IF(G.LT.0D0) CALL ENDJOB(
     >      ' SBR CHAMC :  problem  with  Gradient (=GE+GS-1) < ',I0)
 
C--------- limited to D2G due to IDE = 2
          DG  = (DGE  * QE(1,1)*DSEA + DGS*QS(1,1)*DSSA)  /BR
          D2G = (D2GE * QE(1,2) + D2GS * QS(1,2))  /BR

C--------- Er and derivatives w.r.t R and A
          ER = - G * VN / R
          DRR = - ER / R
          DRA = - DG * RM / R

C--------- Ea and derivatives w.r.t R and A
          RRM = (R-RM) / RM
C          EA = - DG * RRM * ( 1.D0 + ( - 1.5D0 + ( 11.D0 / 6.D0
C     >       - 50.D0/24.D0 * RRM) * RRM) * RRM)
C          DAA = - D2G * (R-RM) * ( 1.D0 + ( - 1.5D0 + ( 11.D0 / 6.D0
C     >       - 50.D0/24.D0 * RRM) * RRM) * RRM)
          DAR = - DG / RM * ( 1.D0 + ( - 3.D0 + ( 5.5D0
     >       - 50.D0/6.D0 * RRM) * RRM) * RRM)
          EA = DRA - R * DAR

C---------- Transform to cartesian co-ordinates
C          Ex, Ey
          E(1,1) = EA
          E(1,2) = ER

          DADX = 1D0 / R

C         ... dBx/dX
C          DE(1,1) = 
C         ... dBx/dY = dBy/dX
C          DE(2,1) = DAR
C         ... dBy/dY
C          DE(2,2) = DRR
  
CC         ... d2Bx/dX2
C          DDE(1,1,1) = ( D2AA2 / R + DRR ) / R       !??? 
CC         ... d2Bx/dXdY = d2By/dX2
C          DDE(2,1,1)= D2ARA * DADX
CC         ... d2Bx/dY2 = d2By/dXdY
C          DDE(2,2,1) = D2RRA * DADX
CC         ... d2By/dY2
C          DDE(2,2,2) = D2RR2
C 
        ENDIF
      ENDIF

      RETURN
      END
