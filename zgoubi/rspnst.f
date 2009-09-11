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
C  Fran�ois M�ot <meot@lpsc.in2p3.fr>
C  Service Acc�lerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      SUBROUTINE RSPNST(MLB,
     >                      PRLB,IA,LBL,NLB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PRLB
      CHARACTER*(*) LBL(MLB)
C     ------------------------------------------
C     Read data for SPNPRNL, SPNSTORE procedures
C     ------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      CHARACTER*80 TA
      COMMON/DONT/ TA(MXL,20)

      INTEGER DEBSTR,FINSTR
      CHARACTER TXT*80, STRA(1)*80

C----- Will print into file TA(NOEL,1), 
C           right after any occurence of element label[s] TA(noel,2)

      READ(NDAT,FMT='(A)') TXT
      TXT = TXT(DEBSTR(TXT):FINSTR(TXT))
      CALL STRGET(TXT,1,
     >                   IDUM,STRA) 
      TA(NOEL,1) = STRA(1)

      IF(MLB .NE. 0) THEN
        ITXT = finstr(txt)
        txt = txt(debstr(txt):itxt)
        LENG = 1+FINSTR(STRA(1))-DEBSTR(STRA(1))
        TXT = TXT(LENG+2:itxt)
        TA(NOEL,2) = TXT
        CALL STRGET(TXT,MLB,
     >                       NLB,LBL)

        PRLB = (NLB .GE. 1) 
     >  .AND. (TA(NOEL,1).NE.'none') .AND. (LBL(1).NE.'none')

        READ(NDAT,*) A(NOEL,1)
C------- Will print every IA turn
        IA=A(NOEL,1)
      ENDIF
      RETURN
      END