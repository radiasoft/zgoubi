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
      SUBROUTINE FITMM(IT,Y,T,Z,P,SAR,DP,TAR,PAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      PARAMETER (MXLOC=100)
      DIMENSION FMI(MXJ,MXLOC,mxt), FMA(MXJ,MXLOC,mxt)
      DIMENSION FLDMI(3,MXLOC,mxt), FLDMA(3,MXLOC,mxt)
      DIMENSION SFLD(3,MXLOC,MXT)

      SAVE FMI, FMA, FLDMI, FLDMA, SFLD

      INCLUDE "MXLD.H"
      DIMENSION IQ(MXL)
      SAVE IQ

      SAVE NBL
      DATA NBL / 0 /

      CALL ZGNOEL(
     >            NOEL)
      NEL = IQ(NOEL)

      FMA(1,NEL,IT) = DMAX1( DP,FMA(1,NEL,IT))
      FMA(2,NEL,IT) = DMAX1(  Y,FMA(2,NEL,IT))
      FMA(3,NEL,IT) = DMAX1(  T,FMA(3,NEL,IT))
      FMA(4,NEL,IT) = DMAX1(  Z,FMA(4,NEL,IT))
      FMA(5,NEL,IT) = DMAX1(  P,FMA(5,NEL,IT))
      FMA(6,NEL,IT) = DMAX1(SAR,FMA(6,NEL,IT))
      FMA(7,NEL,IT) = DMAX1(TAR,FMA(7,NEL,IT))

      FMI(1,NEL,IT) = DMIN1( DP,FMI(1,NEL,IT))
      FMI(2,NEL,IT) = DMIN1(  Y,FMI(2,NEL,IT))
      FMI(3,NEL,IT) = DMIN1(  T,FMI(3,NEL,IT))
      FMI(4,NEL,IT) = DMIN1(  Z,FMI(4,NEL,IT))
      FMI(5,NEL,IT) = DMIN1(  P,FMI(5,NEL,IT))
      FMI(6,NEL,IT) = DMIN1(SAR,FMI(6,NEL,IT))
      FMI(7,NEL,IT) = DMIN1(TAR,FMI(7,NEL,IT))

      FLDMA(1,NEL,IT) = DMAX1(B(1,1),FLDMA(1,NEL,IT))
      FLDMA(2,NEL,IT) = DMAX1(B(1,2),FLDMA(2,NEL,IT))
      FLDMA(3,NEL,IT) = DMAX1(B(1,3),FLDMA(3,NEL,IT))

      FLDMI(1,NEL,IT) = DMIN1(B(1,1),FLDMI(1,NEL,IT))
      FLDMI(2,NEL,IT) = DMIN1(B(1,2),FLDMI(2,NEL,IT))
      FLDMI(3,NEL,IT) = DMIN1(B(1,3),FLDMI(3,NEL,IT))

      SFLD(1,nel,IT)       = SFLD(1,nel,IT) + B(1,1)*PAS
      SFLD(2,nel,IT)       = SFLD(2,nel,IT) + B(1,2)*PAS
      SFLD(3,nel,IT)       = SFLD(3,nel,IT) + B(1,3)*PAS

      RETURN

      ENTRY FITMM1(KT,JI,LI,MIMA,IC2,
     >                               VAL)
C KT=prtcl #, JI=field coordinate, LI=lmnt #
      NLI = IQ(LI)      

      IF    (IC2.GE.1 .AND. IC2.LE.3) THEN
        IF(MIMA.EQ.1) THEN
          VAL = FMI(JI,NLI,KT)
        ELSEIF(MIMA.EQ.2) THEN
          VAL = FMA(JI,NLI,KT)
        ENDIF

      ELSEIF(IC2.GE.6 .AND. IC2.LE.8) THEN
        IF(MIMA.EQ.1) THEN
          VAL = FLDMI(JI,NLI,KT)
        ELSEIF(MIMA.EQ.2) THEN
          VAL = FLDMA(JI,NLI,KT)
        ENDIF

      ELSEIF(IC2.EQ.9) THEN
        
        VAL = SFLD(JI,NLI,KT)
      ELSE
        CALL ENDJOB('SBR FF. NO SUCH CONSTRAINT 3.',IC2)
      ENDIF 

      RETURN

      ENTRY FITMM2(ITI)
      DO  LL = 1, MXLOC
        DO  JJ = 1, MXJ
          FMI(JJ,LL,ITI) = 1.D99
          FMA(JJ,LL,ITI) = -1.D99
        ENDDO
        DO  I = 1, 3
          SFLD(I,LL,ITI) = 0.d0
          FLDMI(i,LL,ITI) = 1.D99
          FLDMA(i,LL,ITI) = -1.D99
        ENDDO
      ENDDO
      RETURN

      ENTRY FITMM4(ic3)
      NBL = NBL + 1
      IQ(ic3) = nbl
      RETURN

      ENTRY FITMM6(I0)
      NBL = I0
      RETURN

      END      
