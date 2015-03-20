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
      SUBROUTINE EVENT(DL,Y,T,Z,P,X,XAR,QBR,SAR,TAR,KEX,IT,
     > AMT,Q,BORO,KART,IMAX,*)
C     > AMT,Q,BORO,KART,IFDES,KGA,KSYN,IMAX,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.AIM_3.H"     ! COMMON/AIM/ ATO,AT,ATOS,RM,XI,XF,EN,EB1,EB2,EG1,EG2
      INCLUDE "MAXTRA.H"
      INCLUDE "C.CHAMBR.H"     ! COMMON/CHAMBR/ LIMIT,IFORM,YLIM2,ZLIM2,SORT(MXT),FMAG,YCH,ZCH
 
      INCLUDE "C.DESIN.H"     ! COMMON/DESIN/ FDES(7,MXT),IFDES,KINFO,IRSAR,IRTET,IRPHI,NDES
C     >,AMS,AMP,AM3,TDVM,TETPHI(2,MXT)
      INCLUDE "C.GASC.H"     ! COMMON/GASC/ AI, DEN, KGA
      INCLUDE "C.SYNRA.H"     ! COMMON/SYNRA/ KSYN

C----- in-flight decay
      IF(IFDES .EQ. 1)
     >  CALL MCDES(DL,KEX,Y,T,Z,P,X,QBR,SAR,TAR,IT,AMT,Q,BORO,XAR,KART)
C     >  CALL MCDES(DL,KEX,Y,T,Z,P,RZ,BR,SAR,TAR,IT,AMT,QT,BORO,XAR,KART)

C----- gas-scattering
      IF(KGA .EQ. 1) THEN
        QBR0 = QBR
        CALL GASCAT(DL,QBR0,IT,
     >                         QBR,*99)
      ENDIF

C------- Walls (chamber) 
      IF(LIMIT .EQ. 1) CALL CHMBRE(IT,Y,Z,SAR,
     >                                        KEX,*99)

C----- synchrotron radiation
      IF(KSYN .GE. 1 ) CALL RAYSYN(DL,IT,IMAX,Q)

      RETURN
 99   RETURN 1
      END
