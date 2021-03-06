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
C  USA
C  -------
      FUNCTION BINARI(NOMFIC,
     >                       IB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) NOMFIC

      LOGICAL BINARI
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG

      INTEGER DEBSTR,FINSTR

      BINARI=.FALSE.

      IDEB = DEBSTR(NOMFIC)
      IFIN = FINSTR(NOMFIC)

      IF( IFIN .LE. IDEB ) RETURN

      IB = IFIN-1
 1    CONTINUE
        IF ( IB .LE. IDEB ) GOTO 2
        IF ( NOMFIC(IB-1:IB-1) .EQ. ']'
     >  .OR. NOMFIC(IB-1:IB-1) .EQ. '/' ) GOTO 2
        IB=IB-1
      GOTO 1

 2    CONTINUE
      BINARI=NOMFIC(IB:IB+1) .EQ. 'B_' .OR. NOMFIC(IB:IB+1) .EQ. 'b_'
      RETURN
      END
