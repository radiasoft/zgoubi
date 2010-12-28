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
      SUBROUTINE SPINR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------
C     Spin rotator. 
C     ------------------
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'       ! MXL, MXD
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXTRA.H"     ! MXT
      INCLUDE "MAXCOO.H"     ! MXJ
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
      COMMON/SPIN/ KSPN,KSO,SI(4,MXT),SF(4,MXT)

      double precision OT(4), VEC(3)

      IOPT = NINT(A(NOEL,1))
      ANG = A(NOEL,10)

      IF(IOPT .NE. 1)THEN
        IF(NRES.GT.0) THEN 
          WRITE(NRES,FMT='(35X,''+++++ OFF +++++'')') 
        ENDIF
        RETURN
      ENDIF

         DO I=1,IMAX
           IF(IEX(I) .ge. 0) then
             S1 = SF(1,i)
             S2 = SF(2,i)
             S3 = SF(3,i)
              dtr=PI/180.d0
              ANGLE=int(ANG)
              phi=(ANG-int(ANG))*100.d0
              VEC(1)=cos(dtr*phi)
              VEC(2)=sin(dtr*phi)
              VEC(3)=0.d0
              Call SPINRO(ANGLE,VEC,OT)
            SF(1,i) = S1*(OT(1)**2+OT(2)**2-OT(3)**2-OT(4)**2)+
     +              S2*2.d0*(OT(2)*OT(3)-OT(1)*OT(4))+
     +              S3*2.d0*(OT(2)*OT(4)+OT(1)*OT(3))
            SF(2,i) = S2*(OT(1)**2-OT(2)**2+OT(3)**2-OT(4)**2)+
     +              S1*2.d0*(OT(2)*OT(3)+OT(1)*OT(4))+
     +              S3*2.d0*(OT(3)*OT(4)-OT(1)*OT(2))
            SF(3,i) = S3*(OT(1)**2-OT(2)**2-OT(3)**2+OT(4)**2)+
     +              S2*2.d0*(OT(3)*OT(4)+OT(1)*OT(2))+
     +              S1*2.d0*(OT(2)*OT(4)-OT(1)*OT(3))
           endif
         enddo

      IF(NRES.GT.0) WRITE(NRES,109) ANGLE, PHI
 109  FORMAT(/,30X,'Spin rotator. Angle = ',F12.5,' deg., ',
     >       /,30X,'              axes at ',F12.5,' deg. ',/)


      RETURN
      END
      
      
