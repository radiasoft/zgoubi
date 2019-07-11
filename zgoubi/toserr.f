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
      SUBROUTINE TOSERR(NOEL,IRR,MXTA,BM,
     >KPOL,TYPERR,TYPAR,TYPDIS,ERRCEN,ERRSIG,ERRCUT,
     >                                     DB,DPOS,TILT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BM(*)
      PARAMETER (JPOL=6)
      DIMENSION KPOL(MXTA,JPOL)
      CHARACTER(2) TYPERR(MXTA,JPOL)
      CHARACTER(1) TYPAR(MXTA,JPOL),TYPDIS(MXTA,JPOL)
      DIMENSION ERRCEN(MXTA,JPOL),ERRSIG(MXTA,JPOL),ERRCUT(MXTA,JPOL)

      INCLUDE 'MXLD.H'
      DIMENSION DB(MXL,JPOL),DPOS(MXL,JPOL,3),TILT(MXL,JPOL,3)

      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST_2.H"     ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      INCLUDE "C.CONST2.H"     ! COMMON/CONST2/ ZERO, UN

      DO I = 1, JPOL

        IF(KPOL(IRR,I) .EQ. 1) THEN

          IF    (TYPDIS(IRR,I).EQ.'G') THEN
            SM = EXP(-(ERRCUT(IRR,I)*ERRSIG(IRR,I))**2/2.D0)
            DERR = (2.D0*RNDM() -1.D0)*SM
          ELSEIF(TYPDIS(IRR,I).EQ.'U') THEN
            XXX = RNDM()
            DERR = ERRSIG(IRR,I)* 2.D0*(XXX-0.5D0)
C            DERR = ERRSIG(IRR,I)* 2.D0*(rndm()-0.5D0)
C                WRITE(88,*) ' TOSERR RNDM ',XXX
          ENDIF

          IF    (TYPERR(IRR,I)(1:1).EQ.'B') THEN
            IF    (TYPAR(IRR,I).EQ.'A') THEN
C              ABSOLUTE ERROR
              DB(NOEL,I) = ERRCEN(IRR,I) + DERR
C              WRITE(NRES,FMT='(10X,A,I,A,2(A,I0),A,1P,3E14.6)')
C     >        'Element #',NOEL,', error added - ',
C     >        ' Multipole order is ',I0,' ;  type of error is ',IRR,
C     >        ' ;  error center/delta/center+delta : ',
C     >        ERRCEN(IRR,I), DERR, DB(NOEL,I)
            ELSEIF(TYPAR(IRR,I).EQ.'R') THEN
C              RELATIVE ERROR
              DB(NOEL,I) = ERRCEN(IRR,I) + DERR*BM(I)
            ENDIF
          ELSEIF(TYPERR(IRR,I)(2:2).EQ.'S') THEN
              DPOS(NOEL,I,3) = 0.D0
          ELSEIF(TYPERR(IRR,I)(2:2).EQ.'R') THEN
              TILT(NOEL,I,3) = 0.D0
          ELSE
            WRITE(NRES,FMT='(//,5X,A)') '****** Found TYPERR = '//TYPERR
            CALL ENDJOB('SBR TOSERR. NO SUCH OPTION FOR TYPERR',-99)
          ENDIF

C        ELSE
C          CALL ENDJOB('SBR TOSERR. NO SUCH OPTION KPOL=',KPOL(IRR,I))

        ENDIF

      ENDDO

      RETURN
      END
