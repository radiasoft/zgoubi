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
C  François Méot <meot@lpsc.in2p3.fr>
C  Service Accélerateurs
C  LPSC Grenoble
C  53 Avenue des Martyrs
C  38026 Grenoble Cedex
C  France
      FUNCTION BONUL(XL,PAS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL BONUL
C     -----------------------------------------------------
C     If  KFLD=0 (no field) , ELEMENT <=> DRIFT
C     -----------------------------------------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "MAXCOO.H"
      INCLUDE "MAXTRA.H"
      LOGICAL AMQLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),IMAX,IEX(MXT),IREP(MXT),AMQLU
      LOGICAL ZSYM
      COMMON/OPTION/ KFLD,MG,LC,ML,ZSYM
 
      BONUL = KFLD .EQ. 0 

      IF(BONUL) THEN
        XXL=XL
        IF(PAS .LT. 0.D0) XXL = -XL
        IF(NRES.GT.0) WRITE(NRES,109) XXL
 109    FORMAT(/,25X,' ++++++++  ZERO  FIELD  ++++++++++' ,/,15X,
     >  'Element  is  equivalent  to  drift  of  length ',F12.5,' cm')
C 109    FORMAT(/,25X,' ++++++++  Champ NUL  ++++++++++' ,/,15X,
C     >  'ELEMENT EQUIVALENT A UN ESPACE  LIBRE  DE',F12.5,' CM')
        CALL ESL(0,XXL,1,IMAX)
      ENDIF
 
      RETURN
      END
