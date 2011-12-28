        subroutine bb
        implicit double precision (a-h,o-z)
      COMMON/CDF/ IES,LF,LST,NDAT,NRES,NPLT,NFAI,NMAP,NSPN,NLOG
      COMMON/CONST/ CL9,CL ,PI,RAD,DEG,QE ,AMPROT, CM2M
      INCLUDE 'MXLD.H'
      COMMON/DON/ A(MXL,MXD),IQ(MXL),IP(MXL),NB,NOEL
      INCLUDE "MAXTRA.H"
      INCLUDE "MAXCOO.H"
      LOGICAL AMQLU(5),PABSLU
      COMMON/FAISC/ F(MXJ,MXT),AMQ(5,MXT),DP0(MXT),IMAX,IEX(MXT),
     $     IREP(MXT),AMQLU,PABSLU
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
        dimension  spincoord(3)
        dimension  sigma(6)
        double precision  coef

        parameter (mxpt=10000)
        dimension  Ptsl(9,mxpt)
        double precision linearmap(6,6) 

        double precision  ga, bet

        save Npt,coef,sigma,Ptsl,ga
        save linearmap

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

          mass=amq(1,1)
          charge=amq(2,1)
          P0 = BORO*CL9*Q
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
        
        open(unit=2,file='zgoubbi.out')
        do i = 1, Npt
           write(2,fmt='(2i6,6(1x,e15.6))') npt,i,(Ptsl(j,i),j=1,6)
        enddo
        close(2)


        open(9,file="ptcl.out",status="unknown",position="append",
     >                                        form="formatted")
   
        if(iturn.eq.1) then
          write(9,*)Npt
        endif
        write(9,*)iturn
        do i =1, Npt
          snorm = sqrt(Ptsl(7,i)**2+Ptsl(8,i)**2+Ptsl(9,i)**2)
          write(9,222) Ptsl(1:9,i),snorm,i,iturn
        enddo

        close(9)
222     format(10(1x,e17.9),2(1x,i6))
        return
        end
