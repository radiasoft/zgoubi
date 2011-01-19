C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  Fran�ois M�ot
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
C  Fran�ois M�ot <fmeot@bnl.gov>
C  Brookhaven National Laboratory               �s
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE IMPVAR(IUNIT,NI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXV=40) 
      COMMON/CONTR/VAT(MXV),XI(MXV)
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      COMMON /VAR/ X(3*MXV),P(MXV)
      COMMON/VARY/NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,7)
      INCLUDE 'MXFS.H'
      PARAMETER (LBLSIZ=8)
      PARAMETER (KSIZ=10)
      CHARACTER FAM*(KSIZ),LBF*(LBLSIZ),KLEY*(KSIZ),LABEL*(LBLSIZ)
      COMMON/SCALT/ FAM(MXF),LBF(MXF,2),KLEY,LABEL(MXL,2)
      CHARACTER*(KSIZ) KLE
      CHARACTER*(LBLSIZ) LBL1, LBL2
      LOGICAL EMPTY

        IF(IUNIT .EQ. 7) WRITE(IUNIT,100)
100     FORMAT('1')
        WRITE(IUNIT,200) NI
200     FORMAT(/,' STATUS OF VARIABLES  (Iteration #',I6,')')
        WRITE(IUNIT,300)
300     FORMAT(
     >  ' LMNT  VAR  PARAM  MINIMUM     INITIAL         FINAL        ',
C------   IR(I)  I   IS(I)   X(K)        XI(I)          X(I)  
     >  ' MAXIMUM      STEP     NAME       LBL1     LBL2' )
C----      X(J)        P(I)
 
      DO 1 I=1,NV
        K=I+NV
        J=K+NV
        CALL ZGKLE(
     >             KLE,IQ(IR(I)))
        LBL1 = LABEL(IR(I),1)
        LBL2 = LABEL(IR(I),2)
        IF(EMPTY(LBL1)) LBL1 = '*'
        IF(EMPTY(LBL2)) LBL2 = '*'
        WRITE(IUNIT,400) IR(I),I,IS(I),X(K),XI(I),
     >  A(IR(I),IS(I)),X(J),P(I),KLE,LBL1,LBL2
 400    FORMAT(1P, 
     >  2X,I3,3X,I2,4X,I3,2(2X,G10.3),2X,G17.10,2(1X,G10.3),3(1X,A))
        IF(XCOU(I).NE.0.D0) THEN
          KL=XCOU(I)
          KP=NINT((1D3*XCOU(I)-1D3*KL))
          SGN = KL/ABS(KL)
          ISGN = NINT(SGN)
          WRITE(IUNIT,400) ISGN*KL,I,ISGN*KP,
     >            A(IR(I),IS(I)),XI(I),SGN*X(I),X(J),P(I)
        ENDIF
 1    CONTINUE
      CALL FLUSH2(IUNIT,.FALSE.)
      RETURN
      END
