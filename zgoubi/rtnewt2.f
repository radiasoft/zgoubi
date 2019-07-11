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
C  USA
C  -------
      SUBROUTINE RTNEWT2(XB,YB,D,E,X1,X2,XACC,TTARF,TZ,FA,Typ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (JMAX=10000)
      PARAMETER(lunY=12,lunW=13)
      INCLUDE "C.SPIRALE.H"     ! COMMON/spiral_ent/UMEG,ASP0,ASP1,ASP2,ASP3
      INCLUDE "C.SPIRALX.H"     ! COMMON/spiral_ext/UMEGs,ASPS0,ASPS1,ASPS2,ASPS3
      INCLUDE "C.RADIALS.H"     ! COMMON/radial_sec/aen,ben,cen,aex,bex,cex

      pi = 4.d0 * atan (1.d0)


       open(unit=lunW,file='plot_spiral.H')

      if (Typ .EQ. 1.0) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCC     Malek plot spiral boundaries   CCCCCCCCCCCCCCCCCCCCCCCC
       if (i .LT. 2) then
            rmin=50.0
            rmax=1000.0
            theta_min=-20.0/180.0*pi        ! -2*pi/12.0
            theta_max=20.0/180.0*pi        ! 2*pi/12.0
            deltat=0.0005
            Nzi=int((theta_max-theta_min)/deltat)

            do iz=1,Nzi+1
               theta=theta_min+(iz-1)*deltat
               r0=rmin+(iz-1)*deltar
               thet=theta
               call RSOLVE(thet,D,rmin,rc,FA)

               if (FA==1.0) then
                  a0=UMEG
               else
                  a0=UMEGS
               endif
               r=rc
               xpro=r*cos(thet+a0)
               ypro=r*sin(thet+a0)
               Write(lunW,*) xpro,ypro
            enddo
        endif
         i=i+1


       elseif (Typ .EQ. 0.0) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC    Magnet plot radial boundaries  CCCCCCCCCCCCCCCCCCCCCCCC

         if (i .LT. 3) THEN
            theta_min=-30.0/180.0*pi        ! -2*pi/12.0
            theta_max=30.0/180.0*pi        ! 2*pi/12.0
            deltat=0.00005
            Nzi=int((theta_max-theta_min)/deltat)

            if (FA==1.0) then
              a=aen
              b=ben
              c=cen
            else
              a=aex
              b=bex
              c=cex
            endif

            do iz=1,Nzi+1
               theta=theta_min+(iz-1)*deltat
               thet=theta

               r=-c/(a*cos(thet)+b*sin(thet))

               xpro=r*cos(thet)
               ypro=r*sin(thet)
!               Write(lunY,*) xpro,ypro
            enddo
          endif
          i=i+1
       endif


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    solve spiral sector  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        if (Typ .EQ. 1.0) THEN

      TZ=X2
C      sol0: i have to figure out how to give a good first approximation to r_sol
      sol0=SQRT(XB**2+YB**2)
      DO 11 J=1,JMAX
         TEMP=TZ
         CALL RSOLVE(TZ,D,sol0,RS,FA)
           IF (FA==1.0) THEN
              xi=ASP0+ASP1*RS+ASP2*RS**2+ASP3*RS**3
           ELSE
              xi=ASPS0+ASPS1*RS+ASPS2*RS**2+ASPS3*RS**3
           ENDIF
         E=1.D0/tan(xi)




CCCCCCCCCCCCCCCCCCCCCC

         CALL FUNCD2(TZ,XB,YB,D,E,F,DF,TTARF)
         DX=F/DF
C        write(*,*) DX
C     DEB RAJOUT
         DXLIM=1.D0
         IF (ABS(DX).GT.DXLIM) THEN
            DX=0.1D0*(X1+X2)
         ENDIF
C     FIN RAJOUT
         TZ=TZ-DX
         IF((X1-TZ)*(TZ-X2).LT.0.D0) THEN
           TZ=TEMP
           RETURN
         ENDIF

         IF(ABS(DX).LT.XACC) RETURN

 11      CONTINUE
C      PAUSE 'rtnewt2 exceeding maximum iterations'

         RETURN

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    solve radial sector  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCC   D = RM,  (XB,YB) = point courant,

       elseif (Typ .EQ. 0.0) THEN


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      endif

      close(lunW)

      close(lunY)

      RETURN
      END
