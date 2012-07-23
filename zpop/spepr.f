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
C  Brookhaven National Laboratory  
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  -------
      SUBROUTINE SPEPR(NLOG,KPR,NT,NPTS,YM,YPM,YNU,PMAX,NC0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YM(*), YPM(*), YNU(*), PMAX(*), NC0(*)

      CHARACTER*12 HVL(3)
      CHARACTER*2 YC(3), YPC(3)
      CHARACTER TIT(3)*4, REP
      CHARACTER TXTP*5,TXTL*14

      DATA YC / 'Y', 'Z', 'X' /
      DATA YPC / 'Y''', 'Z''', 'D' /
      DATA HVL / 'Horizontal', 'Vertical', 'Longitudinal' /
      DATA TIT/'NuY=','NuZ=','NuX='/
 
      IF(NT.EQ.-1) THEN
        WRITE(TXTP,FMT='(A5)') '* all'
      ELSE
        WRITE(TXTP,FMT='(I5)') NT
      ENDIF
      CALL READC3(KL1,KL2)
      IF(KL1.EQ.-1) THEN
        WRITE(TXTL,FMT='(A5)') '* all'
      ELSE
        WRITE(TXTL,FMT='(I5,A,I5)') KL1,' to ',KL2
      ENDIF
      WRITE(*,101) TXTP,TXTL,NPTS 
 101  FORMAT(/,' Part',A5,'  at Lmnt ',A5,' ; ',I6,' PNTS') 

      DO 1 JNU = 1, 3
        WRITE(*,FMT='(A,''  motion :'')') HVL(JNU)
        WRITE(*,FMT=
     >  '(10X,''Center at '',2A,''  ='',1P,2E12.4,'' (MKSA)'')') 
     >  YC(JNU), YPC(JNU), YM(JNU), YPM(JNU)
        WRITE(*,179) TIT(JNU),YNU(JNU),1.D0-YNU(JNU),PMAX(JNU),NC0(JNU)
 179    FORMAT(1X,A4,1P,'/[1-]',G14.6,'/',G14.6,'   Ampl. =',G12.4,
     >        ',   ',I4,' bins')
 1    CONTINUE

      IF(KPR.EQ.1) THEN   
 20     WRITE(*,*)
        WRITE(*,*) '  PRINT IN zpop.log (Y/N)?' 
        READ(*,FMT='(A1)',ERR=20) REP 
        IF(REP .NE. 'N' .AND. REP .NE. 'n') REP = 'y'
      ELSEIF(KPR.EQ.2) THEN
        REP = 'y'
      ENDIF

      IF(REP.EQ. 'y') THEN

        WRITE(*,*) '  Tunes will be printed in zpop.log'
        WRITE(*,*) '---------------------------------------------------'
        WRITE(*,*) 
        WRITE(NLOG,101) TXTP,TXTL,NPTS 

        DO 10 JNU = 1, 3
          WRITE(NLOG,FMT='(A,''  motion :'')') HVL(JNU)
          WRITE(NLOG,FMT=
     >    '(10X,''Center at '',2A,''  ='',1P,2E12.4,'' (MKSA)'')') 
     >    YC(JNU), YPC(JNU), YM(JNU), YPM(JNU)
          WRITE(NLOG,179) 
     >         TIT(JNU),YNU(JNU),1.D0-YNU(JNU),PMAX(JNU),NC0(JNU)
C----------  This is to flush the write statements...
             CALL FLUSH2(NLOG,.FALSE.)
 10     CONTINUE

      ENDIF

      RETURN
      END
