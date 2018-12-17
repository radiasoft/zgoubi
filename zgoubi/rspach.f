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

      SUBROUTINE RSPACH(MLB,
     >                      TSPCH1,TLAMBDA1,LBL,NLB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TSPCH1
      CHARACTER(*) LBL(*)
C     ---------------------------------
C     READS DATA FOR SPACECHARGE PROCEDURE (SIMILAR TO FAISTORE)
C     ---------------------------------
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      LOGICAL TSPCH
      INCLUDE "C.SPACECHA.H" ! COMMON/SPACECHA/ TLAMBDA,RBEAM(2),XAVE(2),EMITT(2),TAVE,BUNCH_LEN,
C                             >                EMITTZ, BTAG, SCKX, SCKY, TSPCH

      INTEGER DEBSTR,FINSTR
      CHARACTER(80) TXT, STRA(1)

c      save TLAMBDA1
c      save TSPCH1

      READ(NDAT,FMT='(A)') TXT
      CALL STRGET(TXT,1,
     >                  IDUM,STRA)
C File name
      TA(NOEL,1) = STRA(1)

      IF(MLB .EQ. 0) THEN
        NLB = 0
        LBL(1) = ' '
        LBL(2) = ' '
        TSPCH1 = .FALSE.
      ELSE
        ITXT = finstr(txt)
        TXT = TXT(DEBSTR(TXT):ITXT)
        LENG = 1+FINSTR(STRA(1))-DEBSTR(STRA(1))
        TXT = TXT(LENG+2:itxt)
        TA(NOEL,2) = TXT
        CALL STRGET(TXT,MLB,
     >                      NLB,LBL)
        TSPCH1 = ((NLB .GE. 1)
     >  .AND. (TA(NOEL,1).NE.'none') .AND. (LBL(1).NE.'none')
     >  .OR. LBL(1).EQ.'all'
     >  .OR. LBL(1).EQ.'ALL' )

        READ(NDAT,*) A(NOEL,1)
        TLAMBDA1=A(NOEL,1)

      ENDIF

      TSPCH = TSPCH1
      TLAMBDA = TLAMBDA1

      RETURN

c       ENTRY SPA1(
c     >             TLAMBDA1,TSPCH1)
c       TLAMBDA1=TLAMBDA
c       TSPCH1=TSPCH
c       RETURN

c       write(*,*) TSPCH1, TLAMBDA1

      END
