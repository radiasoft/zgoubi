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
      SUBROUTINE ELCYLF(MPOL,EM,QLEM,QLSM,QE,QS,A,R,Z,
     >                                          E,DE,DDE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(MCOEF=6)
      DIMENSION QLEM(MPOL),QLSM(MPOL),QE(MPOL,MCOEF),QS(MPOL,MCOEF)
      DIMENSION EM(*)
      DIMENSION E(5,3),DE(3,3),DDE(3,3,3)

      COMMON/ORDRES/ KORD,IRD,IDS,IDB,IDE,IDZ

      DIMENSION D3EX(3,3,3), D3EY(3,3,3), D3EZ(3,3,3)
      DIMENSION D4EX(3,3,3,3) ,D4EY(3,3,3,3) ,D4EZ(3,3,3,3)

C------ ER0(n) = -d^nV/dR^n at Z=0
      PARAMETER (MDR=5)
      DIMENSION ER0(MDR)
      DATA ER0 /  MDR*0.D0 /

c        write(*,*) ' elcylf ',A,R,Z
      CALL ELCYLM(MPOL,EM,QLEM,QLSM,QE,QS,A,R,Z,
     >                                        E,DE,DDE)
c        write(*,*) ' elcylf E ',e
c      CALL EREXYZ(ER0,Z,IDE,
c     >                       E,DE,DDE)
      CALL DBDXYZ(IDE,DE,DDE,D3EX,D3EY,D3EZ,D4EX,D4EY,D4EZ)

      RETURN
      END
