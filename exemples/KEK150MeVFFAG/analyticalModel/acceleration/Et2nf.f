C Transform the  orbit duration vs Energy dependence,  as obtained from eg prior run of 
C geneMap, into freq as a function of turn number. 
C Input data are read in searchCO.out_COs ; Et2nf.In gives additional working 
C conditions data, output is printed Et2nf.out, for further read by the SCALING procedure 
C when running CAVITE for acceleration in zgoubi. 

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      parameter(nCO2=2002,nCO = nCO2/2)
C tau is turn duration on closed orbit at energy E
      dimension turn(nCO2), freq(nCO2), Ekin(nCO2), tau(nCO2)

      character let*1
      integer*4 nbTrn
C-------------------------- 
C User's working data  
      open(unit=7,file='Et2nf.In')  
      read(7,*) nCell, Vp, phsD, E1, E2
      close(7)

        ak = 7.6
      am = 938.27231e6
      kTrStp = 1
C--------------------------

      pi = 4.d0 * atan(1.d0)
      phs = phsD /180. * pi
      dE = Vp  * sin(phs)
      nbTrn = (E2-E1)/dE
C the prgrm works with MeV
      am = am * 1.d-6
      dE = dE * 1.d-6
      E1 = E1 * 1.d-6
      E2 = E2 * 1.d-6  * 1.2

      write(6,*) '  Vp, phsD, E1, E2*1.2, dE, #turns :', 
     >                        Vp, phsD, E1, E2, dE, nbTrn

      open(unit=7,file='Et2nf.In2')  
      open(unit=8,file='Et2nf.Out')  

          i = 1
          j = nCO2
 1    continue
C read closed orbit coordinates and duration (tau)
        read(7,*,err=598,end=599) x,xp,z,xp,s,d,tco,let,E
          if(j.eq.nCO2/2) stop '  Too many data'
          turn(j) = 1.d0+(E-E1)/dE
          tau(j) = tco * nCell 
          freq(j) = 1.d0/tau(j) 
          Ekin(j) = E 
c          write(88,fmt='(1p,3e14.6,2x,i6,a)') 
c     >      turn(j),freq(j),tau(j),j,'  turn#, freq.(MHz), tau(mu_s), j'
          j = nCO2-i
          i = i+1
        goto 1

      close(7)
      close(8)

c            write(88,*) '% -------------------- '

 598      write(*,*) '  read ended upon error, # data read is ',i-1
          goto 21
 599      write(*,*) '  read ended upon EOF, # data read is ',i-1
          goto 21
 21       continue

c            write(88,*) '% -------------------- '

            n = i-1
            write(*,*) ' # of data to be interpolated :  n= ',n
            if(n.eq.0) stop ' SBR CUBSPL, no data to be interpolated'

c        write(88,*) '% -------------------- '

C Reverses the array so to have turn increasing with increasing index
            if(turn(nCO2).gt.turn(nCO2-1)) then
              do 3 i = 1, n
                turn(i) = turn(nCO2-n+i)        
                freq(i) = freq(nCO2-n+i)        
                tau(i) = tau(nCO2-n+i)        
                Ekin(i) = Ekin(nCO2-n+i)        
 3            continue
            else
              do 5 i = 1, n
                turn(n-i+1) = turn(nCO2-n+i)        
                freq(n-i+1) = freq(nCO2-n+i)        
                tau(n-i+1) = tau(nCO2-n+i)        
                Ekin(n-i+1) = Ekin(nCO2-n+i)        
 5           continue
            endif
              do 31 i = 1, n
                write(88,fmt='(1p,4e16.8,2x,i5,2X,a)') 
     >           turn(i), freq(i), tau(i)/nCell, Ekin(i), i,
     >           'turn, freq, tauCell, Ekin, i  '
 31          continue

        write(*,*) '% ---------------- ',n,nint(turn(1)),int(turn(n))
             pause
C------ In order to compute "one" :
        trn = 1.d0
        yv = CUBSPL(turn,freq,trn,nCO2,n)
        tt1 = CUBSPL(turn,tau,trn,nCO2,n)
        ek1 = CUBSPL(turn,Ekin,trn,nCO2,n)
        p1 = ((ek1+am)**2 - am*am)**(0.5d0)
C----------------------------------
              
C------ This is for checking interpolation
      do 61 i = 1, n
        Ek = Ekin(i)
        tt = CUBSPL(Ekin,tau,Ek,nCO2,n)
        frqncy = CUBSPL(Ekin,freq,Ek,nCO2,n)
                write(88,fmt='(1p,4e16.8,2x,i5,2X,a)') 
     >           turn(i), freq(i), tau(i)/nCell, Ekin(i), i,
     >           'turn, freq, tauCell, Ekin, i   '
        write(89,fmt='(1p,2e16.8,2x,i4,2x,5e16.8,a)') tt/nCell, Ek, i,
     >     1.d0/frqncy/nCell-1.d0/freq(i)/nCell, tt-tau(i), turn(i), 
     >       frqncy, freq(i), 
     >      ' tauCell, Ekin, co#'
 61   continue
c      do 62 i = 1, n 
c        Ek = Ekin(i)
c        tt = CUBSPL(Ekin,tau,Ek,nCO2,n)
c        frqncy = CUBSPL(Ekin,freq,Ek,nCO2,n)
c        write(87,fmt='(1p,3e16.8,2x,i3,2x,a)') 
c     >   tau(i), Ekin(i), tt, i,   ' tau, tt, i'
c 62   continue
C----------------------------------

      ntrn2 = int(turn(n))

      sytdx = 0.d0      
      stt = 0.d0
      do 6 ntrn = 1, ntrn2
        trn = dble(ntrn)
        frq = CUBSPL(turn,freq,trn,nCO2,n)
        tt = CUBSPL(turn,tau,trn,nCO2,n)
        stt = stt + tt 
        sytdx = sytdx + frq * tt 
        if(ntrn.eq.1) sytdx = 0.d0        
        phase2 = 2.d0 * pi * sytdx
        phase = 2.d0 * pi * (trn-1.d0)
        ek = CUBSPL(turn,Ekin,trn,nCO2,n)
        p = ((ek+am)**2 - am*am)**(0.5d0)
C           write(*,*) tt,tt1, p, p1, ek, am
C So to check scaling of time duration : 
        one = tt1/tt * (p/p1)**(-ak/(ak+1.d0))  * (ek+am)/(ek1+am)
        PHI = phase + phs
        PHI = PHI - INT(PHI/(2.D0*PI)) * 2.D0*PI 
        IF    (PHI .GT.  PI) THEN
          PHI =PHI - 2.D0*PI
        ELSEIF(PHI .LT. -PI) THEN
          PHI =PHI + 2.D0*PI
        ENDIF
C Print frequency law trn, frq :
        if(mod(ntrn-1,kTrstp) .eq. 0) 
     >   write(8,fmt='(1p,6e16.8,4x,a)') trn, 1.d0/frq/nCell,
     >            phase, stt, ek,
     >  one, !!phi, !!frq*tt, !!one, !! (phase2-phase)/phase, !!    frq*tt, !!
     > 'turn#, tauCell, phi, oclock, Ekin, one'
 6    continue

      stop
      end

      FUNCTION CUBSPL_lin(xm,ym,xv,nd,n)
C Fake
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      doubleprecision XM(nd),YM(nd),LM(nd)

      i = interv(xm,xv,n)

      a = (xm(i+1) - xv) / (xm(i+1) - xm(i))
      b = 1.d0 - a

      cubspl = a * ym(i) + b * ym(i+1)

C       write(*,*) xm(i),xv,xm(i+1),i

      return
      end
      FUNCTION interv(xm,xv,n)
C Fake
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension XM(*)

      if(xv.le.xm(2)) then 
         interv = 1
      elseif(xv.ge.xm(n-1)) then
         interv = n-1
      else
        do jj = 1, n
          if(xm(jj).ge.xv) then 
            j = jj-1
            goto 77
          endif
        enddo
      endif

      return

 77   continue
      interv = j

      return
      end


***CUBIC SPLINE********************************************************
*                  INTERPOLATION USING CUBIC SPLINES                  *
***********************************************************************
*
      FUNCTION CUBSPL(xm,ym,xv,nd,n)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer nd,n,i,m,j,k,count,sele
      doubleprecision 
     >XM(nd),YM(nd),LM(nd),UM(nd),DM(nd),CM(0:nd),EM(nd),
     >xv,fact,e,ff,g,h
** 
**
      cubspl = 9999.d0

      m = n - 1
      j = m - 1
      k = j - 1
**
*  GENERATION OF TRIDIAGONAL SYSTEM FOR SECOND DERIVATIVE
**
      do 10 i = 1,j
         DM(i) = 2.d0* (XM(i+2) - XM(i))
         CM(i) = 6.d0* (YM(i+2) - YM(i+1)) / (XM(i+2) - XM(i+1)) + 6.d0* 
     $                         (YM(i) - YM(i+1)) / (XM(i+1) - XM(i))
10    continue
      do 20 i = 2,j
         LM(i) = XM(i+1) - XM(i)
20    continue
      do 30 i = 1,k
         UM(i) = XM(i+2) - XM(i+1)
30    continue 
**
*  SOLUTION OF TRIDIAGONAL SYSTEM
**
      CALL TRIDI(LM,DM,UM,CM,nd,j)
** 
**
*  EVALUATION AND PRINTING OF CUBIC SPLINES
**
      CM(0) = 0.d0
      CM(n) = 0.d0
c      write(6,*)'------------------------------------------------------'
c      sele = 1
c40    if (sele.eq.1) then
c         write(6,*) ' Enter the value where interpolation is required:'
c         read(*,*) xv
c         write(6,*) 'The equation for cubic splines are:'
         do 50 i = 1,m
            fact = XM(i+1) - XM(i)
            e = CM(i-1) / (6.*fact)
            ff = CM(i) / (6.*fact)
            g = (YM(i)/fact) - (CM(i-1)*fact/6.)
            h = (YM(i+1)/fact) - (CM(i)*fact/6.)
            EM(i) = e* (XM(i+1) - xv) **3 + ff* (xv - XM(i)) **3 + 
     $              g* (XM(i+1)- xv) + h* (xv - XM(i))
c            write(6,41) 'f',i,'x =',e,'(',XM(i+1),'-x)**3 +',ff,
c     $      '(x-',XM(i),')**3 +',g,'(',XM(i+1),'-x)+',h,'(x-',XM(i),')'
c41          format('0',a1,i1,a3,f7.3,a1,f7.3,a8,f7.3,a3,f7.3,a6,f7.3,a1,
c     $             f7.3,a4,f7.3,a3,f7.3,a1)
50       continue 
**
*  SELECTION OF APPROPRIATE SEGMENT (BASED ON THE VALUE) WHERE
*  INTERPOLATION REQUIRED
**
         count = 1
CC FM Oct 2007   60       if (xv.lt.XM(count+1)) go to 70
60       if (xv.le.XM(count+1)) go to 70
         count = count + 1
         go to 60
70       continue 
** 
**
c         write(6,*)'--------------------------------------------------'
c         write(6,71) xv,EM(count)
c71       format('0',5x,'The interpolated value at',f6.2,' is:',f8.3)
c         write(6,*)'--------------------------------------------------'
c         write(6,*)' More Interpolation ?  < 0- for NO / 1- for YES >:'
c         read(*,*) sele
c         go to 40
        CUBSPL = EM(COUNT)
c      endif
      return
c      stop
      end
*
***TRIDIAGONAL MATRIX******************************************************
*                                                                         *
      subroutine TRIDI(L,D,U,COE,nd,n)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      integer nd,n,m,i
      integer m,i
c      parameter(nd=40)
      doubleprecision L(nd),D(nd),U(nd),COE(nd)
**
      m = n - 1 
      do 10 i = 1,m
         L(i+1) = L(i+1) / D(i)
         D(i+1) = D(i+1) - L(i+1) *U(i)
         COE(i+1) = COE(i+1) - L(i+1) *COE(i)
10    continue
**
*  THE COEFFICIENT VECTOR WILL TRANSFORM TO SOLUTION VECTOR*
**
      COE(n) = COE(n)/D(n)
      do 20 i = m,1,-1
         COE(i) = (COE(i) - U(i) *COE(i+1)) / D(i)
20    continue
      return
      end


