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
C  Brookhaven National Laboratory                                                               és
C  C-AD, Bldg 911
C  Upton, NY, 11973
C  USA
C  -------
      LOGICAL FUNCTION INPECH()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OKECH, OKVAR, OKBIN
      COMMON/ECHL/OKECH, OKVAR, OKBIN
      COMMON/VXPLT/ XMI,XMA,YMI,YMA,KX,KY,IAX,LIS,NB

      CHARACTER REP
      save ityp
      data ityp / 2 /

 30   CONTINUE
      CALL HOMCLR

      WRITE(6,103) XMI,XMA,YMI,YMA, ITYP
 103  FORMAT(/,' * Present  values  of  scale  min - max, Mode :',
     >//,'    Xmin,max, Ymin,max, mode(1-3) : ',//, 
     >1P, 4G17.9,3X,I1//)
      
      WRITE(6,100)
 100  FORMAT('$  * Want  to  change  (Y/N) : ')
      READ(5,200,ERR=30) REP 
 200  FORMAT(A1)

      IF(REP.EQ. 'N' .OR. REP.EQ. 'n') THEN
        IF(XMI.LT. XMA .AND. YMI.LT. YMA) THEN
          CALL TRAXES(XMI,XMA,YMI,YMA,ityp)
          OKECH=.TRUE.
        ELSE
          OKECH = .FALSE.
        ENDIF
      ELSEIF(REP.EQ. 'Y' .OR. REP.EQ. 'y') THEN
 40     CONTINUE
        WRITE(6,*)
        WRITE(6,101) 
 101    FORMAT('  Give  XMI, XMA, YMI, YMA, no grid / grid / no axes (1'
     >     '-3) : ')
        READ(5,*,ERR=40) XMI,XMA,YMI,YMA, ITYP
        if(ityp.lt.1 .or. ityp.gt.3 ) ityp = 2
        IF(XMI.LT. XMA .AND. YMI.LT. YMA) THEN 
C          CALL TXTFBG
          CALL TRAXES(XMI,XMA,YMI,YMA,ITYP) 
          CALL FBGTXT
          OKECH=.TRUE.
        ELSE
          WRITE(6,*) '  Min-Max  error !'
          GOTO 30
        ENDIF        
      ELSE
        GOTO 30
      ENDIF
 
      INPECH=OKECH
      RETURN 
      END
