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
        subroutine bb
C S. White & F. Meot, Jan. 2012
        implicit double precision (a-h,o-z)
      INCLUDE "C.CDF.H"     ! COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      INCLUDE "C.CONST.H"     ! COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      INCLUDE "C.DON.H"     ! COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     >     IREP(MXT),AMQLU,PABSLU
      COMMON/REBELO/ NRBLT,IPASS,KWRT,NNDES,STDVM
      INCLUDE "C.RIGID.H"     ! COMMON/RIGID/ BORO,DPREF,DP,QBR,BRI
      COMMON/UNITS/ UNIT(MXJ)

        parameter (epsilon0 = 8.854187817d-12)
        parameter (r0 = 1.535e-18)
        parameter (Gs = 1.79285d0)

        double precision  charge,mass,energev,intensity
        double precision  alfx,betx,epsnx,epsx
        double precision  alfy,bety,epsny,epsy
        double precision  circ,alfmom,dpp,sigz
        double precision  tunex,tuney,tunez
        double precision  ampx,ampy,ampz
C        dimension  spincoord(3)
        dimension  sigma(6)
        double precision  coef

        parameter (mxpt=10000)
        dimension  Ptsl(9,mxpt)
        double precision linearmap(6,6) 

        double precision  ga, bet

        save Npt,coef,sigma,Ptsl,ga
        save linearmap

        logical ok, idluni

c        open(1,file="zgoubbi.in",status="old")
c        read(1,*)Nturns
c        read(1,*)Npt  
c        if(Npt.gt.mxpt) stop ' Npt > max allowed'
c        read(1,*)charge,mass,energev,intensity
c        read(1,*)alfx,betx,epsnx
c        read(1,*)alfy,bety,epsny
c        read(1,*)sigz,dpp
c        read(1,*)circ,alfmom
c        read(1,*)tunex,tuney,tunez
c        read(1,*)spincoord(1:3)
c        read(1,*)ampx,ampy,ampz
c        close(1)

      if(nint(a(noel,1)).eq.0) then
          write(nres,*) ' BEAM-BEAM is OFF'
         
      else
        if(ipass.eq.1) then

       write(nres,*) '   intensity= ',a(noel,2)
       write(nres,*)  ' alfx= ',a(noel,10)
       write(nres,*)  ' betx= ',a(noel,11)
       write(nres,*)  ' epsnx= ',a(noel,12)
       write(nres,*)  ' alfy= ',a(noel,20)
       write(nres,*)  ' bety= ',a(noel,21)
       write(nres,*)  ' epsny= ',a(noel,22)
       write(nres,*)  ' sigz= ',a(noel,30)
       write(nres,*)  ' dpp= ',a(noel,31)
       write(nres,*)  ' circ= ',a(noel,40)
       write(nres,*)  ' alfmom= ',a(noel,41)
       write(nres,*)  ' tunex= ',a(noel,50)
       write(nres,*)  ' tuney= ',a(noel,51)
       write(nres,*)  ' tunez= ',a(noel,52)
       write(nres,*)  ' ampx= ',a(noel,60)
       write(nres,*)  ' ampy= ',a(noel,61)
       write(nres,*)  ' ampz= ',a(noel,62)

          mass=amq(1,1)
          charge=amq(2,1)
          P0 = BORO*CL9*charge
          energev=.001d0 * sqrt(p0*p0 + mass*mass)
          intensity=a(noel,2)
          alfx=a(noel,10)
          betx=a(noel,11)
          epsnx=a(noel,12)
          alfy=a(noel,20)
          bety=a(noel,21)
          epsny=a(noel,22)
          sigz=a(noel,30)
          dpp=a(noel,31)
          circ=a(noel,40)
          alfmom=a(noel,41)
          tunex=a(noel,50)
          tuney=a(noel,51)
          tunez=a(noel,52)
****        spincoord(1)
          ampx=a(noel,60)
          ampy=a(noel,61)
          ampz=a(noel,62)

          ga = sqrt(energev**2+mass**2)/mass
          bet = sqrt(1.0-1.0/ga/ga)
          epsx = epsnx/ga/bet
          epsy = epsny/ga/bet

          gx = (1.0+alfx*alfx)/betx
          gy = (1.0+alfy*alfy)/bety

          sigma(1) = betx*epsx
          sigma(2) = -alfx*epsx
          sigma(3) = gx*epsx
          sigma(4) = bety*epsy
          sigma(5) = -alfy*epsy
          sigma(6) = gy*epsy

          coef = charge*qe*intensity*(1.0+bet**2) 
     >          /(ga*bet*(bet+bet)*2*pi*epsilon0*mass)

          call bblmap(linearmap,tunex,tuney,tunez,
     >              alfx,betx,alfy,bety,alfmom,ga,bet,cl,circ)

        endif
      endif

      return

      entry bbkick

      call bbkck(Npt,coef/2.0d0,sigma,Ptsl,ga,Gs)
      call bbkLm(Ptsl,Npt,linearmap)
        
      if(ipass.eq.1) then
c        ok = idluni(
c     >              luna)
c        open(unit=luna,file='zgoubbi.out')

        ok = idluni(
     >              lunb)
        open(lunb,file="ptcl.out",status="unknown",position="append",
     >                                        form="formatted")
          write(lunb,*)Npt
      endif

        do i = 1, Npt
c           write(luna,fmt='(2i6,6(1x,e15.6))') npt,i,(Ptsl(j,i),j=1,6)
          f(1,i) = Ptsl(6,i)
          f(2,i) = Ptsl(1,i)*1d2
          f(3,i) = Ptsl(2,i)*1d3
          f(4,i) = Ptsl(3,i)*1d2
          f(5,i) = Ptsl(4,i)*1d3
ccccccc          f(6,i) = Ptsl(5,i)*1d6
        enddo
c        close(luna)


        write(lunb,*) ipass
        do i =1, Npt
          snorm = sqrt(Ptsl(7,i)**2+Ptsl(8,i)**2+Ptsl(9,i)**2)
          write(lunb,222) Ptsl(1:9,i),snorm,i,ipass
        enddo

c        close(lunb)
222     format(10(1x,e17.9),2(1x,i6))
        return
        end
