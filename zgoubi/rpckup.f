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
      SUBROUTINE RPCKUP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      PARAMETER (MXPUD=9,MXPU=1000)
      COMMON/CO/ FPU(MXPUD,MXPU),KCO,NPU,NFPU,IPU
      PARAMETER (MPULAB=5)
      CHARACTER*10 PULAB
      COMMON/COT/ PULAB(MPULAB)

C----- NPU = 0 (OFF) or NPU > 0 (# of distinct lmnt LABEL's)
      READ(NDAT,*) NPU
      IF(NPU .GT. MPULAB) THEN 
        NPU = MPULAB
        WRITE(NRES,FMT='(/,10X'' ** WARNING '', 
     >  ''  Too  many  PUs  =>  number  forced  to  '',I3)') MPULAB 
      ELSEIF(NPU .GT. 0) THEN
        KCO = 1
      ENDIF

C----- LABEL's
      READ(NDAT,*) (PULAB(I),I=1,NPU)

      RETURN
      END
