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
      SUBROUTINE IMPCTR(IUNIT,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXV=60) 
      INCLUDE "C.CONTR.H"     ! COMMON/CONTR/VAT(MXV),XI(MXV)
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      PARAMETER (LBLSIZ=20)
      CHARACTER(LBLSIZ) LABEL
      INCLUDE "C.LABEL.H"     ! COMMON/LABEL/ LABEL(MXL,2)
      INCLUDE "C.FINIT.H"     !COMMON/FINIT/ FINI
      INCLUDE "C.VAR.H"     ! COMMON /VAR/ X(3*MXV),P(MXV)
      INCLUDE "C.VARY.H"  ! COMMON/VARY/ NV,IR(MXV),NC,I1(MXV),I2(MXV),V(MXV),IS(MXV),W(MXV),
                          !     >IC(MXV),IC2(MXV),I3(MXV),XCOU(MXV),CPAR(MXV,27)
      INCLUDE 'MXFS.H'
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) FAM ; CHARACTER(LBLSIZ) LBF
      INCLUDE "C.SCALT.H"     ! COMMON/SCALT/ FAM(MXF),LBF(MXF,MLF)
      CHARACTER(KSIZ) KLE
      CHARACTER(LBLSIZ) LBL1, LBL2
      LOGICAL EMPTY
      LOGICAL NOSYS,NOSYSI,NSYS,NSYSI
      SAVE NOSYS, NSYS
      CHARACTER(108) TXTBUF

      DATA NSYS, NOSYS / .FALSE. , .FALSE. /

      IF(NSYS) THEN
        IF(IUNIT.GT.0) THEN
          WRITE(TXTBUF,FMT='
     >    ('' Fit reached penalty value '',1P,E12.4,
     >     ''            still working on it... '')') FINI
          CALL ARRIER(TXTBUF)
C          CALL FLUSH2(IUNIT,.FALSE.) 
        ENDIF
        RETURN
      ENDIF     

      CALL ZGPNLT( 
     >            PNLTGT)

      IF(IUNIT.GT.0) WRITE(IUNIT,500) PNLTGT
500   FORMAT(' STATUS OF CONSTRAINTS (Target penalty = ',1P,E12.4,')')
      IF(IUNIT.GT.0) WRITE(IUNIT,600)
600   FORMAT(
     >' TYPE  I   J  LMNT#       DESIRED           WEIGHT       ',
     >'  REACHED         KI2       NAME       LBL1     LBL2',
     >'         *  Parameter(s) ')  
C----      X(J)        P(I)
      DO I=1,NC
        XI2=((VAT(I)-V(I))/W(I))**2/F
        NPRM1 = NINT(CPAR(I,1)) + 1
        CALL ZGKLE(IQ(I3(I))
     >                      ,KLE)
        LBL1 = LABEL(I3(I),1)
        LBL2 = LABEL(I3(I),2)
        IF(EMPTY(LBL1)) LBL1 = '*'
        IF(EMPTY(LBL2)) LBL2 = '*'
        IF(IUNIT.GT.0) WRITE(IUNIT,700) 
     >  IC(I),I1(I),I2(I),I3(I),V(I),W(I),VAT(I),XI2, 
     >  KLE,LBL1,LBL2,NINT(CPAR(I,1)),(CPAR(I,JJ),JJ=2,NPRM1)
700     FORMAT(1P,3I4,I6,5X,E14.7,4X,E11.4,3X,E14.7,2X,E11.4,2X,
     >  3(1X,A),' * ',I3,' : ',6(E9.1,'/'))
      ENDDO
      IF(IUNIT.GT.0) THEN
        WRITE(IUNIT,FMT='
     >  ('' Fit reached penalty value '',1P,E12.4)') FINI
        CALL FLUSH2(IUNIT,.FALSE.) 
      ENDIF

      NSYS = NOSYS

      RETURN

      ENTRY IMPCT2(NOSYSI,NSYSI)
      NOSYS = NOSYSI
      NSYS = NSYSI
      RETURN
      END
