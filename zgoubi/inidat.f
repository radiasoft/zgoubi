C  ZGOUBI, a program for computing the trajectories of charged particles
C  in electric and magnetic fields
C  Copyright (C) 1988-2007  François Meot
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
C  François Meot <fmeot@bnl.gov>
C  BNL
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE INIDAT
      USE DYNHC
      use pariz_namelist_interface, only : ID, MXX, MXY, IZ, MMAP
      DOUBLE PRECISION CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M

      INCLUDE "C.CONST.H"  ! COMMON/CONST/ CL9,CL,PI,RAD,DEG,QEL,AMPROT,CM2M
      DATA ISTAT / 0 /

!     Allocate storage for array A accordingly
      IF( .NOT.ALLOCATED( HC ))
     >  ALLOCATE( HC(ID,MXX,MXY,IZ,MMAP), STAT = ISTAT)
      IF (ISTAT .NE. 0)
     >     CALL ENDJOB('SBR INIDAT Not enough memory for Malloc of HC',
     >     -99)


      CL9=CL*1.D-9
      PI = 4.D0 * ATAN(1.D0)
      RAD = PI/180.D0
      DEG = 180.D0/PI
      CM2M = 1.D-2

      CALL INIARR(1)
      CALL KSMAP0

      RETURN

      END
