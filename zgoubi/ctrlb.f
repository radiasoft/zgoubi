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
      SUBROUTINE CTRLB(K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE "C.CNTRLB.H"     ! COMMON/CNTRLB/ XSTP,DIVB,ALAPL(4),RTN(3)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG

      GOTO (1,2) K

 1    CONTINUE
      XSTP = 0
      DIVB = 0D0
      DO 11 I=1,3
        ALAPL(I) = 0.D0
        RTN(I) = 0.D0
 11   CONTINUE
      ALAPL(4)=0.D0

      RETURN

 2    CONTINUE
      IF(XSTP .EQ. 0.D0) RETURN
      DIVB = DIVB/XSTP
      DO 21 I = 1,3
        ALAPL(I) = ALAPL(I)/XSTP
 21     RTN(I) = RTN(I)/XSTP
      ALAPL(4) = ALAPL(4)/XSTP

      IF(NRES.GT.0) WRITE(NRES,103) XSTP,DIVB,
     1(ALAPL(I),RTN(I),I=1,3),ALAPL(4)
 103  FORMAT(//,15X,' CONDITIONS  DE  MAXWELL  (',F10.0,'  PAS )  :'
     1,/,   20X,'   DIV(B)   ',5X,'LAPLACIEN(B)',5X,'ROTATIONNEL(B)'
     1,/,   19X,   1P,G12.4    ,5X,   G12.4   ,6X,  G12.4
     1, 2(/,36X                  ,    G12.4   ,6X,  G12.4)
     1,/,   20X,'   LAPLACIEN SCALAIRE =',G12.4,//)

      RETURN

      END
