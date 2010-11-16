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
C  Brookhaven National Laboratory                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      SUBROUTINE SRPLV(NLOG,KSC,MX,SYDX,FAC,XVA,XDI,YVA,YDI,SY,
     >                                                            YNRM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SYDX(*)
      CHARACTER*(*) XVA,XDI,YVA,YDI,SY

      CHARACTER YDI1*15, SY1*13

      IF(YNRM .EQ. 1.D0) THEN
        YDI1 = YDI
        SY1 = SY
      ELSE      
        YDI1 = '(Normalized)'
        SY = '(Normalized)'
      ENDIF

      CALL TRKVAR(MX,YVA,YDI1,XVA,XDI)
      CALL SRPLI(NLOG,KSC,SYDX,FAC,SY)

      IF(YNRM .EQ. 1.D0) THEN
        WRITE(* ,100)
        WRITE(NLOG,100)
 100    FORMAT(' (normalize by I/e (A/C) to get P(W) in the window.')
      ELSE      
        WRITE(* ,101) YNRM
        WRITE(NLOG,101) YNRM
 101    FORMAT(/,' The spectrum and integrals are normalized by',E16.8)
      ENDIF
      RETURN
      END
