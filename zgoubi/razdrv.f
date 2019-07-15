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
C  USA
C  -------
      SUBROUTINE RAZDRV(IOP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     --------------------------------------------------------
C     INITIALISATIONS DANS CHAMC.
C     INITIALISATIONS LIEES A L'ORDRE DE CALCUL ( 2, 25 OU 4 )
C     DU Champ B(X,Y,Z).
C     --------------------------------------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CHAMP.H"     ! COMMON/CHAMP/ BZ0(5,5), EZ0(5,5)
      INCLUDE "C.CHAVE.H"     ! COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      INCLUDE "C.DDBXYZ.H"     ! COMMON/DDBXYZ/ DB(3,3),DDB(3,3,3)
      INCLUDE "C.D3B.H"     ! COMMON/D3BXYZ/ D3BX(27), D3BY(27), D3BZ(27)
      INCLUDE "C.D4B_2.H"     ! COMMON/D4BXYZ/ D4BX(81) ,D4BY(81) ,D4BZ(81)
      INCLUDE "C.DDEXYZ.H"     ! COMMON/DDEXYZ/ DE(3,3),DDE(3,3,3)
      INCLUDE "C.D3E.H"     ! COMMON/D3EXYZ/ D3EX(27), D3EY(27), D3EZ(27)
      INCLUDE "C.D4EXYZ_2.H"     ! COMMON/D4EXYZ/ D4EX(81), D4EY(81), D4EZ(81)

      GOTO (1,2,1) IOP

 1    CONTINUE

      B(1,:) = 0.D0
      DB(:,:) = 0.D0
      DDB(:,:,:) = 0.D0
      D3BX(:) = 0.D0
      D3BY(:) = 0.D0
      D3BZ(:) = 0.D0
      D4BX(:) = 0.D0
      D4BY(:) = 0.D0
      D4BZ(:) = 0.D0
      BZ0(:,:) = 0.D0

      IF(IOP .EQ. 1) RETURN

 2    CONTINUE

      E(1,:) = 0.D0
      DE(:,:) = 0.D0
      DDE(:,:,:) = 0.D0
      D3EX(:) = 0.D0
      D3EY(:) = 0.D0
      D3EZ(:) = 0.D0
      D4EX(:) = 0.D0
      D4EY(:) = 0.D0
      D4EZ(:) = 0.D0
      EZ0(:,:) = 0.D0

      RETURN
      END
