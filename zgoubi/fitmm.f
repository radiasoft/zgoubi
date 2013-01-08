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
      SUBROUTINE FITMM(Y,T,Z,P,SAR,DP,TAR,PAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CHAVE/ B(5,3),V(5,3),E(5,3)
      INCLUDE "MAXCOO.H"
      PARAMETER (MXLOC=100)
      DIMENSION FMI(MXJ,MXLOC), FMA(MXJ,MXLOC)
      DIMENSION FLDMI(3,MXLOC), FLDMA(3,MXLOC)
      DIMENSION SFLD(3,MXLOC)

      SAVE FMI, FMA, FLDMI, FLDMA, SFLD

      INCLUDE "MXLD.H"
      DIMENSION IQ(MXL)
      SAVE IQ

      SAVE NBL
      DATA NBL / 0 /

      CALL ZGNOEL(
     >            NOEL)
      NEL = IQ(NOEL)

c        WRITE(*,*) IQ
c          READ(*,*)

c      write(88,*) 'fitmm 1 ',NOEL,NEL,Y, FMA(2,NEL), MAX(Y,FMA(2,NEL))
c      write(*,*) 'fitmm 1 ',NOEL,NEL,Y, FMA(2,NEL), MAX(Y,FMA(2,NEL))
      FMA(1,NEL) = MAX(DP,FMA(1,NEL))
      FMA(2,NEL) = MAX(Y,FMA(2,NEL))
      FMA(3,NEL) = MAX(T,FMA(3,NEL))
      FMA(4,NEL) = MAX(Z,FMA(4,NEL))
      FMA(5,NEL) = MAX(P,FMA(5,NEL))
      FMA(6,NEL) = MAX(SAR,FMA(6,NEL))
      FMA(7,NEL) = MAX(TAR,FMA(7,NEL))

      FMI(1,NEL) = MIN(DP,FMI(1,NEL))
      FMI(2,NEL) = MIN(Y,FMI(2,NEL))
      FMI(3,NEL) = MIN(T,FMI(3,NEL))
      FMI(4,NEL) = MIN(Z,FMI(4,NEL))
      FMI(5,NEL) = MIN(P,FMI(5,NEL))
      FMI(6,NEL) = MIN(SAR,FMI(6,NEL))
      FMI(7,NEL) = MIN(TAR,FMI(7,NEL))

      FLDMA(1,NEL) = MAX(B(1,1),FLDMA(1,NEL))
      FLDMA(2,NEL) = MAX(B(1,2),FLDMA(2,NEL))
      FLDMA(3,NEL) = MAX(B(1,3),FLDMA(3,NEL))

      FLDMI(1,NEL) = MIN(B(1,1),FLDMI(1,NEL))
      FLDMI(2,NEL) = MIN(B(1,2),FLDMI(2,NEL))
      FLDMI(3,NEL) = MIN(B(1,3),FLDMI(3,NEL))

      SFLD(1,nel)       = SFLD(1,nel) + B(1,1)*PAS
      SFLD(2,nel)       = SFLD(2,nel) + B(1,2)*PAS
      SFLD(3,nel)       = SFLD(3,nel) + B(1,3)*PAS

c      write(88,fmt='(a,2(1x,i3),3(1x,1p,e12.4))') 'fitmm sfld ',NOEL,NEL
c     > ,sfld(1,nel),sfld(2,nel),sfld(3,nel)
          
      RETURN

      ENTRY FITMM1(JI,LI,MIMA,IC2,
     >                            VAL)
      NLI = IQ(LI)      
      IF    (IC2.GE.1 .AND. IC2.LE.3) THEN
        IF(MIMA.EQ.1) THEN
          VAL = FMI(JI,NLI)
        ELSEIF(MIMA.EQ.2) THEN
          VAL = FMA(JI,NLI)
        ENDIF
      ELSEIF(IC2.GE.6 .AND. IC2.LE.8) THEN
        IF(MIMA.EQ.1) THEN
          VAL = FLDMI(JI,NLI)
        ELSEIF(MIMA.EQ.2) THEN
          VAL = FLDMA(JI,NLI)
        ENDIF
      ELSEIF(IC2.EQ.9) THEN
        VAL = SFLD(JI,NLI)
      ELSE
        CALL ENDJOB('SBR FF. NO SUCH CONSTRAINT 3.',IC2)
      ENDIF 

c      write(*,fmt='(a,3i4,e14.6)') ' FITMM IC2,JI,LI,VAL ',IC2,JI,LI,VAL
c         read(*,*)

      RETURN

      ENTRY FITMM2
c        write(*,*) ' fitmm2 ',mxloc, mxj
c           read(*,*)
      DO  LL = 1, MXLOC
        DO  JJ = 1, MXJ
          FMI(JJ,LL) = 1.D99
          FMA(JJ,LL) = -1.D99
c        write(*,*) ' fitmm2 ',FMI(JJ,LL),FMa(JJ,LL) 
        ENDDO
        DO  I = 1, 3
          SFLD(I,LL) = 0.d0
          FLDMI(i,LL) = 1.D99
          FLDMA(i,LL) = -1.D99
        ENDDO
      ENDDO
      RETURN

      ENTRY FITMM3(ic3)
      NBL = NBL + 1
      IQ(ic3) = nbl
c        write(*,*) ' fitmm nbl, ic3, iq(ic3) ',nbl, ic3, iq(ic3)
c         read(*,*)
      RETURN

      END      
