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
C  Brookhaven National Laboratory                és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE DATE2(DMY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER * 9   DMY

C This works :   ------------------------------- 
C      CHARACTER*3 MM(12)
C      CHARACTER*8 DD
C      CHARACTER*10 TT
C      CHARACTER*5 ZZ
C      INTEGER VV(8)
C      DATA MM/'Jan','Feb','Mar','Apr','May','Jun',
C     $     'Jul','Aug','Sep','Oct','Nov','Dec'/
C      LONG = LEN(DMY)
C      IF(LONG.GE.9) THEN
C         CALL DATE_AND_TIME(DD,TT,ZZ,VV)
C         WRITE(DMY,'(I2.2,''-'',A3,''-'',I2.2)')
C     $        VV(3),MM(VV(2)),MOD(VV(1),100)
C      ENDIF
C      IF(LONG.GT.9) DMY=DMY(1:9)//' '      
C-------------------------------------------------

C This uses built-in idate, hence transportable
      INTEGER TODAY(3)
      call idate(today)  ! today(1)=day, today(2)=month, today(3)=year
      write(dmy,fmt='(2(i2.2,a1),i2.2)') 
     > today(1),'/',today(2),'/',today(3) !-2000
      RETURN
      END
