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
      FUNCTION APDPO(A,IRA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=500)
      DIMENSION FP(N)
      DATA EPS / 1.D-6 /
      SAVE EPS
      IA = INT(A)
      IF(4*IA.GE.N) THEN
        WRITE(6,FMT='(/,5X,
     >  '' ** WARNING : You need larger array FP for calculation of'',
     >  '' poisson law'',/,10X,''Suggestion : decrease step size'')')
        STOP
      ENDIF
CC----- Calculate normalized POISSON law p(k)=exp(-A)*A^k/k!
      NN=N
      CALL APDPO1(A,EPS,NN,FP,I1)
C         NN now contains the index k of the last non-zero value
      FPMAX= -1.D-10
      JJ=-1
      INTA=1+IA
 11   JJ=-JJ
      INTA = INTA+JJ
      J=INTA
      IF(J.GT.NN) GOTO 11
 10   CONTINUE
        FPMAX = FP(J)
        JMAX=J
        J=J+JJ
        IF(J.EQ.0) GOTO 12
        IF(J.GT.NN) GOTO 11
        IF(FP(J).GE.FPMAX) THEN
          GOTO 10
        ELSEIF(J.EQ.INTA+1) THEN
          GOTO 11
        ENDIF
        IF(J.EQ.2+INTA) GOTO 11
 12     CONTINUE
C      SUM=0.D0
      DO 15 J=1,NN
        FPJ=FP(J)
        FP(J)=FPJ/FPMAX
        IF(FPJ.LT.0.D0) THEN
C--------- 'cause sometimes FP(NN) slightly < 0 !
           WRITE(6,fmt='(/,A,I6,A)') 
     >     ' ** WARNING / poisson law : p(k)<0 at k=',J,
     >     '    action undertaken : extrapolated from p(k-2), p(k-1)'
           FP(J)=FP(J-1)/(FP(J-2)/FP(J-1))
        ENDIF
 15     CONTINUE
C-------- FP(J) now contains the normalised POISSON law
C            FP_k, k=0,nn-1 = FP(J=1,NN)
C----- Number of photons emitted by each particle, 
C                 "acceptance-rejection" method. 
 17     R1=RNDM(IRA)       
          K= NN*R1
          J=K+1
          R2=RNDM(IRA)       
        IF(R2.GT.FP(J)) GOTO 17
       APDPO=K
       RETURN
       END
