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
      SUBROUTINE RSRLOS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
C      PARAMETER (LNTA=132) ; CHARACTER(LNTA) TA
C      PARAMETER (MXTA=45)
      INCLUDE "C.DONT.H"     ! COMMON/DONT/ TA(MXL,MXTA)
      PARAMETER (MSTR=10)
      CHARACTER(80) STRA(MSTR)
      CHARACTER(132) TXT132
      INTEGER DEBSTR, FINSTR
      LOGICAL STRCON
      PARAMETER (KSIZ=10)
      CHARACTER(KSIZ) KLE

      LINE = 1
      READ(NDAT,*,ERR=90,END=90) A(NOEL,1)

      LINE = LINE + 1
      READ(NDAT,FMT='(A80)',ERR=90,END=90) TXT132
      IF(STRCON(TXT132,'!',
     >                     IS)) TXT132 = TXT132(1:IS-1)
      TA(NOEL,1)=TXT132(DEBSTR(TXT132):FINSTR(TXT132))
      CALL STRGET(TA(NOEL,1),10,
     >                         NSTR,STRA)
      IF(NSTR .GT. MSTR) GOTO 90
      TA(NOEL,1)=' '
      TA(NOEL,2)=' '
      IF(NSTR.GE.1) TA(NOEL,1)=STRA(1)
      IF(NSTR.GE.2) TA(NOEL,2)=STRA(2)
      IF(NSTR.GE.3) THEN
C Get the list of elements to be subjected to scaling
        TA(NOEL,3) = STRA(3)(debstr(STRA(3)):finstr(STRA(3)))
        DO I = 4, NSTR
          TA(NOEL,3)=TA(NOEL,3)(debstr(TA(NOEL,3)):finstr(TA(NOEL,3)))
     >    //' '//STRA(I)(debstr(STRA(I)):finstr(STRA(I)))
        ENDDO
      ELSE
        TA(NOEL,3)=' ' 
      ENDIF

      LINE = LINE + 1
      READ(NDAT,*,ERR=90,END=90) A(NOEL,10), A(NOEL,11)

      RETURN

 90   CONTINUE
      CALL ZGKLEY( 
     >            KLE)
      CALL ENDJOB('*** Pgm rsrlos, keyword '//KLE//' : '// 
     >'input data error, at line #',LINE)
      RETURN
      END
