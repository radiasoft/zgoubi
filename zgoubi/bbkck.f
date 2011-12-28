       subroutine bbkck(Npt,coef,sigma,Ptsl,ga,Gs)
       implicit double precision (a-h,o-z)
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
           c = 0
           sa = 1 - cos(phi)
           sb = sin(phi)
 
           Ptsl(6,i)=Ptsl(6,i)-bbgx*sigma(2)+ 
     >      bbgy*sigma(5)
           Ptsl(6,i)=Ptsl(6,i)-(bbfx*(Ptsl(2,i)-bbfx*0.5)+
     >      bbfy*(Ptsl(4,i)-bbfy*0.5))*0.5
           Ptsl(2,i)=Ptsl(2,i)-bbfx
           Ptsl(4,i)=Ptsl(4,i)-bbfy
           call bbrsp(Ptsl,i,a,b,c,sa,sb,npt)

         enddo
      return
      end 
