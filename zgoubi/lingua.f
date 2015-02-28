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
C  Brookhaven National Laboratory                    és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE LINGUA(LNG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      CHARACTER   TITRE*80 
      INTEGER DEBSTR
 
      READ (NDAT,103) TITRE
 103  FORMAT(A)
      IDEB=DEBSTR(TITRE)
      IF( TITRE(IDEB:IDEB) .EQ. '*') THEN
        IF( TITRE(IDEB+2:IDEB+2) .EQ. '*') THEN
          IF    ( TITRE(IDEB+1:IDEB+1) .EQ. 'F') THEN
C----------- FRANCAIS
            LNG=1
          ELSEIF( TITRE(IDEB+1:IDEB+1) .EQ. 'E') THEN
C----------- ENGLISH
            LNG=2
          ELSEIF( TITRE(IDEB+1:IDEB+1) .EQ. 'J') THEN
C----------- NIPPON
            LNG=3
          ELSEIF( TITRE(IDEB+1:IDEB+1) .EQ. 'V') THEN
C----------- VIET
            LNG=4
          ENDIF
        ENDIF
      ELSE
        LNG=1
      ENDIF
      REWIND(NDAT)
      RETURN
      END
