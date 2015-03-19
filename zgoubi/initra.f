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
      SUBROUTINE INITRA(I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ***COORDONNEES DE LA TRAJECTOIRE I ,INITIALISATION
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      INCLUDE "C.FAISC.H"     ! COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
C     $     IREP(MXT),AMQLU,PABSLU
      INCLUDE "C.PTICUL.H"     ! COMMON/PTICUL/ AM,Q,G,TO
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/TRAJ/ Y,T,Z,P,X,SAR,TAR,KEX,IT,AMT,QT


C///////////////////////
c      call ZGNOEL(
c     >             NOEL)
c       if(noel.eq.245) then
c           write(*,*) 'initra IN : ',noel,i
c           write(*,*) ' QBR,DP ',QBR,DP
c           write(*,*) '  F(1,I)  ',F(1,I)
c         endif

      IT = I
      Y=F(2,I)
      T=F(3,I)*0.001D0
      Z=F(4,I)
      P=F(5,I)*0.001D0
      DP=F(1,I)
      QT = AMQ(2,I)
      QBR = Q*BORO*DP

      BRI = QT/QBR
      KEX=IEX(I)
      SAR= F(6,I)
      AMT = AMQ(1,I)
C----- AMQ(2,I) = Q/QE
      TAR = F(7,I)   *1.D5

cC///////////////////////
c       if(noel.eq.245) then
c           write(*,*) 'initra OUT : ',noel,i
c           write(*,*) ' QBR,DP ',QBR,DP
c           write(*,*) ' BRI  ',bri
c         endif

      RETURN
      END
