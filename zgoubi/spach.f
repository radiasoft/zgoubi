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
      SUBROUTINE SPACH(KSPCH,LBLSC,NLBSC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) LBLSC(*)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL

      KSPCH = NINT(A(NOEL,1))

      IF(NRES.GT.0) THEN 
        IF    (KSPCH .NE. 1) THEN
          WRITE(NRES,107)
 107      FORMAT(/,15X,' KSPCH .ne. 1 :  Space charge is  off. ',/)
        ELSE
          WRITE(NRES,FMT=
     >    '(/,15X,''  KSPCH = 1 :  Space charge switched on. '',/)')

          WRITE(NRES,FMT='(/,A,I0,A,/)') 
     >    'Space charge will be applied at the ',NLBSC,
     >    ' following labels : '

          DO ILBL = 1, NLBSC
            WRITE(NRES,FMT='(5X,A)') LBLSC(ILBL)
          ENDDO

        ENDIF
      ENDIF
 
      RETURN

      ENTRY SCKICK

      RETURN
      END
