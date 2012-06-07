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
      SUBROUTINE FITMM(Y,T,Z,P,SAR,DP,TAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "MAXCOO.H"
      INCLUDE "MXLD.H"
      DIMENSION FMI(MXJ,MXL), FMA(MXJ,MXL)

      DIMENSION FMIO(MXJ,MXL), FMAO(MXJ,MXL)
      SAVE FMIO, FMAO

      CALL ZGNOEL(
     >            NOEL)

      FMA(1,NOEL) = MAX(DP,FMA(1,NOEL))
      FMA(2,NOEL) = MAX(Y,FMA(2,NOEL))
      FMA(3,NOEL) = MAX(T,FMA(3,NOEL))
      FMA(4,NOEL) = MAX(Z,FMA(4,NOEL))
      FMA(5,NOEL) = MAX(P,FMA(5,NOEL))
      FMA(6,NOEL) = MAX(SAR,FMA(6,NOEL))
      FMA(7,NOEL) = MAX(TAR,FMA(7,NOEL))
      FMI(1,NOEL) = MIN(DP,FMI(1,NOEL))
      FMI(2,NOEL) = MIN(Y,FMI(2,NOEL))
      FMI(3,NOEL) = MIN(T,FMI(3,NOEL))
      FMI(4,NOEL) = MIN(Z,FMI(4,NOEL))
      FMI(5,NOEL) = MIN(P,FMI(5,NOEL))
      FMI(6,NOEL) = MIN(SAR,FMI(6,NOEL))
      FMI(7,NOEL) = MIN(TAR,FMI(7,NOEL))
          
      DO JJ = 1, MXJ
          FMIO(JJ,NOEL) = FMI(JJ,NOEL) 
          FMAO(JJ,NOEL) = FMA(JJ,NOEL)
C         write(*,*) jj,noel,FMIO(JJ,NOEL),FMaO(JJ,NOEL)
      ENDDO

C        write(*,*) ' fitmm y,t,noel ',y,t,noel
C        write(*,*) ' fitmm y,t,noel ',FMAO(2,NOEL),FMAO(3,NOEL),noel

      RETURN

      ENTRY FITMM1(JI,LI,MIMA,
     >                        VAL)
      IF(MIMA.EQ.1) THEN
        VAL = FMIO(JI,LI)
      ELSEIF(MIMA.EQ.2) THEN
        VAL = FMAO(JI,LI)
      ENDIF
C        write(*,*) ' fitmm ji,li,mima,val ',ji,li,mima,val
C         stop
      RETURN

      ENTRY FITMM2
      DO  LL = 1, MXL
        DO  JJ = 1, MXJ
          FMI(JJ,LL) = 1.D99
          FMA(JJ,LL) = -1.D99
        ENDDO
      ENDDO
      RETURN

      END      
