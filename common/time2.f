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
C  Brookhaven National Laboratory                �s
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE TIME2(HMS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER * 9  HMS

C This works :   ------------------------------- 
C      CHARACTER*8 DD
C      CHARACTER*10 TT
C      CHARACTER*5 ZZ
C      INTEGER VV(8)
C      LONG = LEN(HMS)
C      IF(LONG.GE.8) THEN
C         CALL DATE_AND_TIME(DD,TT,ZZ,VV)
C         WRITE(HMS,'(I2.2,'':'',I2.2,'':'',I2.2)') VV(5),VV(6),VV(7)
C      ENDIF
C      IF(LONG.GT.8) HMS=HMS(1:8)//' '
c-----------------------------------------------

C This uses built-in idate, hence transportable
      INTEGER NOW(3)
      call itime(now)  ! now()=hour,   now()=minute,   now()=secomd
      write(hms,fmt='(2(i2.2,a1),i2.2)') now(1),':',now(2),':',now(3)


      RETURN
      END
