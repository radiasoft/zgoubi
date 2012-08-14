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
       subroutine bbkck(Npt,coef,sigma,Ptsl,ga,Gs)
       implicit double precision (a-h,o-z)
C S. White & F. Meot, Jan. 2012
       double precision  coef,sepx,sepy,bbfx,bbfy,bbgx,bbgy
       dimension  sigma (6)
       parameter (mxpt=10000)
       dimension  Ptsl(9,mxpt)
       double precision  pxdz,pydz,phi
       double precision  a,b,c,sa,sb
       double precision  ga, Gs

         do i=1,npt

           sepx=Ptsl(1,i)
           sepy=Ptsl(3,i)

           call bb4G
     >        (sepx,sepy,sigma(1),sigma(4),bbfx,bbfy,bbgx,bbgy)

           bbfx=-coef*bbfx
           bbfy=-coef*bbfy
           bbgx=-coef*bbgx
           bbgy=-coef*bbgy

           pxdz = -bbfx*((1+ga*Gs)+(ga*Gs+ga/(1+ga)))/2.0
           pydz = bbfy*((1+ga*Gs)+(ga*Gs+ga/(1+ga)))/2.0
           phi = sqrt(pxdz**2+pydz**2)

           a = pxdz/phi
           b = pydz/phi
           c = 0.d0
           sa = 1 - cos(phi)
           sb = sin(phi)
 
           Ptsl(6,i)=Ptsl(6,i)-bbgx*sigma(2)+ 
     >      bbgy*sigma(5)
           Ptsl(6,i)=Ptsl(6,i)-(bbfx*(Ptsl(2,i)-bbfx*0.5)+
     >      bbfy*(Ptsl(4,i)-bbfy*0.5))*0.5
           Ptsl(2,i)=Ptsl(2,i)-bbfx
           Ptsl(4,i)=Ptsl(4,i)-bbfy
           call bbrsp(Ptsl,i,a,b,c,sa,sb)
C           call bbrsp(Ptsl,i,a,b,c,sa,sb,npt)

         enddo
      return
      end 
