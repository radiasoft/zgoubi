        subroutine bb4G(sepx,sepy,sigxx,sigyy,bbfx,bbfy,
     >              bbgx,bbgy)
        implicit double precision (a-h,o-z)
c        implicit none
        double precision pieni,x,const,xxyy,expfac,bbfx,bbfy,
     >                  sepx,sepy,sigxx,sigyy,bbgx,bbgy,comfac,comfac2

        pieni = 1.d-38
        
        x=sepx**2+sepy**2
        xxyy=sigxx+sigyy
        const=0.0d0
        if(abs(xxyy).gt.pieni) const=x/xxyy
        expfac=exp(-const)
        bbfx=0.0d0
        bbfy=0.0d0
        bbgx=0.0d0
        bbgy=0.0d0
        if(abs(x).gt.pieni) then          
          bbfx=2.0d0*sepx*(1.0d0-expfac)/x
          bbfy=2.0d0*sepy*(1.0d0-expfac)/x
          comfac=-sepx*bbfx+sepy*bbfy
          comfac2=(abs(sigxx)+abs(sigyy))**2
          bbgx=(comfac+4.0d0*sepx**2*const/x*expfac)/(2.d0*x)
          bbgy=(-comfac+4.0d0*sepy**2*const/x*expfac)/(2.d0*x)
        endif
        return
        end 
