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
      SUBROUTINE SPEGR(NT,JNU,YNU,PMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YNU(*), PMAX(*)

      CHARACTER TIT(3)*5, TXT*80
      PARAMETER (PTXT0=21.D0)
      SAVE PTXT
      CHARACTER TXTP*5,TXTL*14
      DATA TIT/'NuY =','NuZ =','NUX ='/
      DATA PTXT / PTXT0 /

      CALL LOGO

      WRITE(TXT,100) 'TUNES'
 100  FORMAT(A5)
      CALL TRTXT(120.D0,245.D0,TXT,0) 
      WRITE(TXT,179) YNU(JNU),'  amplitude = ',PMAX(JNU)
 179    FORMAT(F10.6,A14,G12.4 )
      TXT = TIT(JNU)//TXT(1:75)
      CALL TRTXT(10.D0,PTXT,TXT,0)
      PTXT = PTXT - 7.D0
      IF(PTXT .LT. 7.D0) PTXT = PTXT0
      IF(NT.EQ.-1) THEN
        WRITE(TXTP,FMT='(A5)') '* all'
      ELSE
        WRITE(TXTP,FMT='(I5)') NT
      ENDIF
      CALL READC3(KL1,KL2)
      IF(KL1.EQ.-1) THEN
        WRITE(TXTL,FMT='(A5)') '* all'
      ELSE
        WRITE(TXTL,FMT='(I5,A5,I5)') KL1,' to ',KL2
      ENDIF
      WRITE(TXT,101) TXTP,TXTL
 101  FORMAT(' Part',A,'  at Lmnt ',A) 
      CALL TRTXT(10.D0,0.001D0,TXT,0)

      RETURN
      END
