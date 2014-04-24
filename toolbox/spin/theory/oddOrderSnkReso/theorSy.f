      implicit double precision (a-h,o-z)
      
      pi = 4.d0 * atan(1.d0)
      pi2 = pi/2.d0
      am = .93827203e9
      eps = 6.8d-2   !   3.88d6*12.1e-10       ! resonance strength = (J_n**2/epsilon)*epsilon
      eps2 = eps * eps
      Vp = 3.d5
      phis = 20.d0 *pi/180.d0
      G = 1.7928474
      Qy = 29.6725d0
      GgR = 411.d0-Qy
      Gg0 =  Ggr - 6.1d0      !!  378.267d0
      dlta = GgR - Gg0
      ggscal = 1.d0
      u = pi * eps * .5d0      
      dW = Vp *sin(phis)
      dGg = ggscal * G/am * dW
      write(*,*) ' dGg = ', dGg
      Gg = Gg0
 1    continue
        dlta2 = dlta * dlta
        almbd = sqrt(dlta2 + eps2)
        v = pi2 * almbd 
        sinc = sin(v)/v
        b = u * sinc
        b2 = b*b
        a2 = 1.d0 - b2
        sy = 1.d0 -8.d0 *a2 * b2
        Gg = Gg + dGg
        write(88,*) Gg, sy
        if(Gg.gt.391.d0) goto 10
        dlta = dlta - dGg
      goto 1 

 10   continue
      write(*,*) ' Job completed !'
      stop
      end 

      
      
      
      





