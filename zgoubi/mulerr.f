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
      SUBROUTINE MULERR(NOEL,IRR,MXTA,BM,
     >KPOL,TYPERR,TYPAR,TYPDIS,ERRCEN,ERRSIG,ERRCUT,
     >                                     DB,DPOS,TILT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BM(*)
C      PARAMETER (MXTA=45)
C      PARAMETER (MXERR=MXTA)
      PARAMETER (MPOL=10)
      DIMENSION KPOL(MXTA,MPOL)
      CHARACTER(2) TYPERR(MXTA,MPOL)
      CHARACTER(1) TYPAR(MXTA,MPOL),TYPDIS(MXTA,MPOL)
      DIMENSION ERRCEN(MXTA,MPOL),ERRSIG(MXTA,MPOL),ERRCUT(MXTA,MPOL)

      INCLUDE 'MXLD.H'
      DIMENSION DB(MXL,MPOL),DPOS(MXL,MPOL,3),TILT(MXL,MPOL,3)

      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN

      DO I = 1, MPOL
        IF(KPOL(IRR,I) .EQ.1) THEN
          IF    (TYPDIS(IRR,I).EQ.'G') THEN
            SM = EXP(-(ERRCUT(IRR,I)*ERRSIG(IRR,I))**2/2.D0)
            DERR = (2.D0*RNDM() -1.D0)*SM
          ELSEIF(TYPDIS(IRR,I).EQ.'U') THEN
            XXX = RNDM()
            DERR = ERRSIG(IRR,I)* 2.D0*(XXX-0.5D0)
C            DERR = ERRSIG(IRR,I)* 2.D0*(rndm()-0.5D0)
C                WRITE(88,*) ' MULERR RNDM ',XXX
          ENDIF
          IF    (TYPERR(IRR,I)(1:1).EQ.'B') THEN
            IF    (TYPAR(IRR,I).EQ.'A') THEN
C              ABSOLUTE ERROR
              DB(NOEL,I) = ERRCEN(IRR,I) + DERR
            ELSEIF(TYPAR(IRR,I).EQ.'R') THEN
C              RELATIVE ERROR
              DB(NOEL,I) = ERRCEN(IRR,I) + DERR*BM(I)
            ENDIF
          ELSEIF(TYPERR(IRR,I)(2:2).EQ.'S') THEN
            IF    (TYPERR(IRR,I)(1:1).EQ.'X') THEN
              DPOS(NOEL,I,1) = ERRCEN(IRR,I) + DERR
            ELSEIF(TYPERR(IRR,I)(1:1).EQ.'Y') THEN
              DPOS(NOEL,I,2) = ERRCEN(IRR,I) + DERR
            ELSEIF(TYPERR(IRR,I)(1:1).EQ.'Z') THEN
              DPOS(NOEL,I,3) = ERRCEN(IRR,I) + DERR
            ENDIF
C            DPOS(NOEL,I,3) = 0.D0
          ELSEIF(TYPERR(IRR,I)(2:2).EQ.'R') THEN
            IF    (TYPERR(IRR,I)(1:1).EQ.'X') THEN
              TILT(NOEL,I,1) = ERRCEN(IRR,I) + DERR
            ELSEIF(TYPERR(IRR,I)(1:1).EQ.'Y') THEN
              TILT(NOEL,I,2) = ERRCEN(IRR,I) + DERR
            ELSEIF(TYPERR(IRR,I)(1:1).EQ.'Z') THEN
              TILT(NOEL,I,3) = ERRCEN(IRR,I) + DERR
            ENDIF
C              TILT(NOEL,I,3) = 0.D0
          ELSE
            CALL ENDJOB('Sbr mulerr. No such option for TYPERR',-99)
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END
