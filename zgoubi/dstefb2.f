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
      FUNCTION DSTEFB2(XB,YB,D,AMIN,AMAX,XACC,TTARF,YN,FA,Typ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/spiral_ent/UMEG,ASP0,ASP1,ASP2,ASP3   
      COMMON/spiral_ext/UMEGs,ASPS0,ASPS1,ASPS2,ASPS3
      COMMON/radial_sec/aen,ben,cen,aex,bex,cex 

      parameter(lunW=12)

C FM Jan 2015
      DSTEFB2 = 999.D0 

cc      open(unit=lunW,file='Distance.H')
CCC later discussion on spiral/radial sector      IF (ASP0 .NE. 0.) THEN
      CALL RTNEWT2(XB,YB,D,E,AMIN,AMAX,XACC,TTARF,TZ,FA,Typ)

      IF (Typ==1.0) THEN
       XN = D*EXP(E*TZ)*COS(TZ+TTARF)
       YN = D*EXP(E*TZ)*SIN(TZ+TTARF)

      DSTEFB2=SQRT( (XB - XN)**2 + (YB - YN)**2 ) 

       
      ELSEIF (Typ==0.0) THEN
         IF (FA==1.0) THEN
            a=aen
            b=ben
            c=cen
            
         ELSE
            a=aex
            b=bex
            c=cex
         ENDIF   


         DSTEFB2=abs(a*XB+b*YB+c)/sqrt(a**2+b**2)
         YN=(a**2*YB-a*b*XB-b*c)/(a**2+b**2)



!       XN = -c/(a*cos(TZ+TTARF)+b*sin(TZ+TTARF))*COS(TZ+TTARF)
!       YN = -c/(a*cos(TZ+TTARF)+b*sin(TZ+TTARF))*SIN(TZ+TTARF) 

      ENDIF
 
!!        write(lunW,*) XB,YB,DSTEFB2

 
  
       RETURN

!!       close(lunW)

      END
