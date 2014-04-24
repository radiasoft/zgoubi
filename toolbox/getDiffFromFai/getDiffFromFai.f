      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXJ=7)
      PARAMETER (MXS=4)
      DIMENSION FO(MXJ),F(MXJ),SI(MXS),SF(MXS)

      logical idluni      

      CHARACTER*1 LET
      CHARACTER*8 LBL1, LBL2
      CHARACTER KLEY*10, txt132*132

          IF (IDLUNI(LR)) THEN
           OPEN(UNIT=LR,FILE='zgoubi.fai',ERR=699)
          ELSE
            GOTO 698
          ENDIF

          IF (IDLUNI(ltmp)) THEN
           OPEN(UNIT=ltmp,FILE='temp.fai',ERR=699)
          ELSE
            GOTO 698
          ENDIF

          IF (IDLUNI(lout)) THEN
           OPEN(UNIT=lout,FILE='getDiffFromFai.out',ERR=699)
          ELSE
            GOTO 698
          ENDIF

C Look for number of praticles in zgoubi.fai --------------------------
C      Read 4-line header
          do i = 1, 4
            READ(LR,fmt='(a)',ERR=99,END=10) txt132
            write(*,*) txt132(1:20)
          enddo
          
          nt = 0
          mt = 0
          it1 = 0
          iq0 = -999
          f1mi = 1.d20
          dpp1 =9999.
 2        CONTINUE
            READ(LR,110,ERR=99,END=10)
     >      KEX,(FO(J),J=1,7),
     >      (F(J),J=1,5),s,time, 
     >      (SI(J),J=1,4),(SF(J),J=1,4),
     >      ENEKI, ENERG, 
     >      IT, IREP, SORT, AMQ1,AMQ2,AMQ3,AMQ4,AMQ5, RET, DPR,  PS,
     >      BORO, IPASS,NOEL, 
     >      TX1,KLEY,TX1,TX1,LBL1,TX1,TX1,LBL2,TX1,TX1,LET,TX1
            INCLUDE "FRMFAI.H"
            if(it.gt.it1) then
              mt = mt+1
              it1 = it
              if(kex .gt.0) then
                nt = nt+1
                if    (abs(f(1)) .lt. 1.d-15) then
                  iq0 = it
                  mt0 = mt
                else
                  if(abs(f(1)) .lt. f1mi) then 
                    f1mi = abs(f(1))
                    dpp1 = f(1)
                    iq1 = it
                    mt1 = mt
                  endif
                endif
              endif
            else
              goto 3
            endif
          goto 2

 3        continue
          write(6,*) '   ' 
          write(6,*) ' -----------------------------------------------' 
          write(6,*) ' There are ',mt,' particles in zgoubi.fai'
          write(6,*) '    including ',nt,' with iex > 0.' 
          write(6,*) '   ' 
          write(6,*) ' Particle number ',mt0,' (it=',iq0,') has dp/p ~0'
          write(6,*) ' Particle number ',mt1,' (it=',iq1,') has dp/p ='
     >       ,dpp1,'  and will be used for computation of dispersion' 
          write(6,*) ' -----------------------------------------------' 
          write(6,*) '   ' 
          if(mt0 .gt.mt1) then
            ita = mt1
            itb = mt0
          else
            itb = mt1
            ita = mt0
          endif

          if(ita.eq.itb) stop ' ita=itb => Cannot compute !!!!!!!!'
          write(6,*) ' ita, itb = ',ita,itb 
          write(6,*) '   ' 
C----------------------------------------------------------------------------


C Read coordinates of those 2 particles used for dispersion calculation --------
       rewind(lr)
C      Read 4-line header
          do i = 1, 4
            READ(LR,fmt='(a)',ERR=99,END=10) txt132
          enddo
          
      sxco = 0.d0
      syco = 0.d0
      sDx = 0.d0
      sDy = 0.d0
      noc = 0
 1    CONTINUE
            if(ita.gt.1) then
              do ii = 1, ita-1 
                READ(LR,110,ERR=99,END=10)              
                
              enddo
            endif
            READ(LR,110,ERR=99,END=10)
     >      KEX,(FO(J),J=1,7),
     >      (F(J),J=1,5),s,time, 
     >      (SI(J),J=1,4),(SF(J),J=1,4),
     >      ENEKI, ENERG, 
     >      IT, IREP, SORT, AMQ1,AMQ2,AMQ3,AMQ4,AMQ5, RET, DPR,  PS,
     >      BORO, IPASS,NOEL, 
     >      TX1,KLEY,TX1,TX1,LBL1,TX1,TX1,LBL2,TX1,TX1,LET,TX1
C            INCLUDE "FRMFAI.H" 
C                itu1 = it       
           if(itb.gt.ita+1) then
              do ii = 1, itb-(ita+1) 
                READ(LR,110,ERR=99,END=10)              
              enddo
            endif
            READ(LR,110,ERR=99,END=10)
     >      KEX,(FO(J),J=1,7),
     >      dpp,x,xp,y,yp,(F(J),J=6,7), 
     >      (SI(J),J=1,4),(SF(J),J=1,4),
     >      ENEKI, ENERG, 
     >      ITc, IREP, SORT, AMQ1,AMQ2,AMQ3,AMQ4,AMQ5, RET, DPR,  PS,
     >      BORO, IPASS,NOEL, 
     >      TX1,KLEY,TX1,TX1,LBL1,TX1,TX1,LBL2,TX1,TX1,LET,TX1
C                itu2 = itc       
            
            if(itb+1.le.mt) then
              do ii = itb+1, mt 
                READ(LR,110,ERR=99,END=10)              
              enddo
            endif

C            write(*,*) '  itu1, itu2 : ', itu1, itu2
            if(ita.eq.mt0) then
              xco = f(2)
              yco = f(4)
              sign = 1.d0
            else
              xco = x
              yco = y
              dpp = f(1)
              sign = -1.d0
            endif

            dltx  = x -f(2)
            dltxp = xp-f(3)
            dlty  = y -f(4)
            dltyp = yy-f(5)

            Dx = dltx / dpp * sign
            Dy = dlty / dpp * sign

            write(lout,fmt='(1p,6e14.6,2i4,a)') 
     >      s,xco,yco,Dx,Dy,dpp,it,itc,
     >      ' s, xco, yco, Dx, Dy, dp/p, it, itc'

            sxco = sxco + xco     
            syco = syco + yco     
            sDx = sDx + Dx     
            sDy = sDy + Dy     
            noc = noc + 1

          goto 1

 10   continue

      sxco = sxco /dble(noc)
      syco = syco /dble(noc)
      sDx = sDx /dble(noc)
      sDy = sDy /dble(noc)

      rewind(lout)
      sx2 = 0.d0
      sy2 = 0.d0
      sDx2 = 0.d0
      sDy2 = 0.d0
      noc2 = 0
 4    continue
        read(lout,*,err=49,end=40) s,xco,yco,Dx,Dy
        sx2 = sx2 + xco * xco
        sy2 = sy2 + yco * yco
        sDx2 = sDx2 + Dx * Dx
        sDy2 = sDy2 + Dy * Dy
        noc2 = noc2 + 1
      goto 4

 40   continue
      sx2 = sx2 / noc2
      sy2 = sy2 / noc2
      sDx2 = sDx2 / noc2
      sDy2 = sDy2 / noc2

      sigx = sqrt(sx2 - sxco*sxco)      
      sigy = sqrt(sy2 - syco*syco)      
      sigDx = sqrt(sDx2 - sDx*sDx)      
      sigDy = sqrt(sDy2 - sDy*sDy)      

      backspace(lout)
      write(lout,fmt='(a)') 
     > ' averages : x, y, Dx, Dy, sigx, sigy, sigDx, sigDy,'//
     >' ; check noc.eq.noc2 !'
      write(lout,fmt='(a,1p,8(1x,e14.6),2(1x,i4))') 
     >'# ',sxco,syco,sDx,sDy,sigx,sigy,sigDx,sigDy,noc,noc2

      stop 'GetDiffFromFai completed !'

  49  CONTINUE
      stop 'Error during read in getDiffFromFai'
  99  CONTINUE
      stop 'Error during read in zgoubi.fai'
 699  CONTINUE
      stop 'Could not open file'
 698  CONTINUE
      stop 'No idle unit '
      end

      FUNCTION IDLUNI(LN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IDLUNI

      LOGICAL OPN

      I = 20
 1    CONTINUE
        INQUIRE(UNIT=I,ERR=99,IOSTAT=IOS,OPENED=OPN)
        I = I+1
        IF(I .EQ. 100) GOTO 99
        IF(OPN) GOTO 1
        IF(IOS .GT. 0) GOTO 1
      
      LN = I-1
      IDLUNI = .TRUE.
      RETURN

 99   CONTINUE
      LN = 0
      IDLUNI = .FALSE.
      RETURN
      END
