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
      SUBROUTINE ITER(A,B,C,DS,COSTA,KEX,IT,*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST2/ ZERO, UN
      COMMON/DEPL/ XF(3),DXF(3),DQBRO,DTAR
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      PARAMETER (EPS=1.D-6,ITMAX=1000)
      PARAMETER (EPS2=1.D-10)

      dm=1.d30
CALCUL INTERSECTION DE LA TRAJECTOIRE AVEC DROITE AX+BY+C=0
      DO 1 I=1,ITMAX
        D=A*XF(1)+B*XF(2)+C
        absd = abs(d)
        IF(ABSD .LE. EPS)  then
C Etienne F. 25 Oct. 05
          IF(ABSD .LE. EPS2)  then
            return
          elseif(d.gt.dm.or.d.eq.0.d0 ) then 
            return
          endif
        endif
        dm=absd
        DS=DS-D/COSTA
        IF(DS*DS .GT. 1.D10) CALL KSTOP(5,IT,KEX,*99)
        CALL DEPLA(DS)
 1    CONTINUE

      CALL KSTOP(5,IT,KEX,*99) 

 99   CONTINUE
      RETURN 1
      END
