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
      SUBROUTINE DIAGQ(nlog,lis,M,N,IP,QXMIN,QXMAX,QYMIN,QYMAX)
C      implicit double precision (a-h,o-z)
      double precision x0
      parameter (x0=0,i0=0)

      ii = abs(m) + abs(n)
 
C*****TRACE DE LA LIGNE M*NUX + N*NUZ = IP
 
      IPLUME=4
C      write(*,*) '*****TRACE DE ',N,'*NUZ=',IP
      IF(M.NE.0) GO TO 1
      Y=IP
      Y=Y/N
      IF ((Y.LT.QYMIN).OR.(Y.GT.QYMAX)) GO TO 10
      CALL VECTPH(QXMIN,Y,4)
      IF(LIS .EQ. 2) write(nlog,*) '  '
      IF(LIS .EQ. 2) CALL IMPV(NLOG,NOC,dble(qxmin),dble(Y),x0,x0,ii)
      CALL VECTPH(QXMAX,Y,2)
      IF(LIS .EQ. 2) CALL IMPV(NLOG,NOC,dble(qxmax),dble(Y),x0,x0,ii)

      GO TO 10
    1 CONTINUE
C      write(*,*) '*****TRACE DE ',M,'*NUX=',IP
      IF(N.NE.0) GO TO 2
      X=IP
      X=X/M
      IF ((X.LT.QXMIN).OR.(X.GT.QXMAX)) GO TO 10
      CALL VECTPH(X,QYMIN,4)
      IF(LIS .EQ. 2) write(nlog,*) '  '
      IF(LIS .EQ. 2) CALL IMPV(NLOG,NOC,dble(x),dble(qymin),x0,x0,ii)
      CALL VECTPH(X,QYMAX,2)
      IF(LIS .EQ. 2) CALL IMPV(NLOG,NOC,dble(x),dble(qymax),x0,x0,ii)
      GO TO 10
C      write(*,*) '*****TRACE DE ',M,'*NUX+',N,'*NUZ=',IP
    2 CONTINUE
C                                     INTERS AVEC X=QXMIN
      Y=(IP-M*QXMIN)/N
      IF ((Y.GE.QYMIN).AND.(Y.LE.QYMAX)) THEN
         CALL VECTPH(QXMIN,Y,IPLUME)
         IF(LIS .EQ. 2 .and. iplume.eq.4) write(nlog,*) '  '
         IF(LIS .EQ. 2) CALL IMPV(NLOG,NOC,dble(qxmin),dble(Y),x0,x0,ii)
         IPLUME=2
      ENDIF
C                                     INTERS AVEC X=QXMAX
      Y=(IP-M*QXMAX)/N
      IF ((Y.GE.QYMIN).AND.(Y.LE.QYMAX)) THEN
         CALL VECTPH(QXMAX,Y,IPLUME)
         IF(LIS .EQ. 2 .and. iplume.eq.4) write(nlog,*) '  '
         IF(LIS .EQ. 2) CALL IMPV(NLOG,NOC,dble(qxmax),dble(Y),x0,x0,ii)
         IPLUME=2
      ENDIF
C                                     INTERS AVEC Y=QYMIN
      X=(IP-N*QYMIN)/M
      IF ((X.GE.QXMIN).AND.(X.LE.QXMAX)) THEN
         CALL VECTPH(X,QYMIN,IPLUME)
         IF(LIS .EQ. 2 .and. iplume.eq.4) write(nlog,*) '  '
         IF(LIS .EQ. 2) CALL IMPV(NLOG,NOC,dble(x),dble(qymin),x0,x0,ii)
         IPLUME=2
      ENDIF
C                                     INTERS AVEC Y=QYMAX
      X=(IP-N*QYMAX)/M
      IF ((X.GE.QXMIN).AND.(X.LE.QXMAX)) THEN
         CALL VECTPH(X,QYMAX,IPLUME)
         IF(LIS .EQ. 2 .and. iplume.eq.4) write(nlog,*) '  '
         IF(LIS .EQ. 2) CALL IMPV(NLOG,NOC,dble(x),dble(qymax),x0,x0,ii)
         IPLUME=2
      ENDIF
C
   10 CONTINUE

      call fbgtxt
      write(*,*) ' ',M,'* QY + ',N,'* QZ = ',IP
      IF(LIS .EQ. 2) write(nlog,*) ' ',M,'* QY + ',N,'* QZ = ',IP

      RETURN
      END
