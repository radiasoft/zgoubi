c      SUBROUTINE SPEANA(YM,BORNE,NC0,
c     >                               YNU,SPEC,PMAX)
      DIMENSION YM(3), BORNE(2)
      parameter(mxpts=100000)
      dimension spin(4,mxpts)
      parameter(mxc=20000)
      dimension SPEC(mxc), q(mxc)
      parameter(mxspec=mxpts/100)
      dimension Qs(mxspec)
      PARAMETER ( PI=3.1415926536 , DEUXPI=2.0*PI )
      logical idluni
      DATA KSX, KSY / 1, 2 /
      data nc0 / 4000 /
      logical okend 
      data okend / .false. /
      data borne / 0.48d0, .60d0 / 
C  will compute tunes over Istp-long intervals, from I1 to I2
      data I1, I2 / 1, mxpts /
      Istp = 1000
      if(Istp.gt.mxpts) Istp=mxpts 

      Nspec = (I2-I1)/(Istp-1)      
      write(*,*) 'Spectrum interval : ',Istp,' turns.'

      if(idluni(lunR)) then
        open(unit=lunR,file='spectrum.In')
      else
        stop ' PGM SPECTRUM : no idle unit'
      endif

      if(idluni(lunW2)) 
     >   open(unit=lunW2,file='spectrum.Out2')

      jspec = 1
      ii = 1
 7    continue

      i = 1
 1    continue
        read(lunR,*,end=3) (spin(j,i),j=1,3)
        spin(4,j) = (spin(1,i)**2 +spin(2,i)**2 +spin(3,i)**2)**0.5d0
C        write(89,*) ' spectrum i :',i
        if(i.ge.Istp) goto 2
        if(i.ge.mxpts) goto 2
        i = i + 1
      goto 1

 3    i = i-1
      okend = .true.
 2    continue
      npt = i

      if(npt.le.2) goto 22
      write(*,*) ' Spectrum # : ',jspec,', points ',(jspec-1)*Istp+1,
     >  ' to ',  jspec*Istp
      
        ANUI = BORNE(1)
        ANUF = BORNE(2)
        DELNU=(ANUF - ANUI) / NC0
        PAS=DEUXPI * DELNU
        VAL=DEUXPI *(ANUI - 0.5d0 * DELNU)
        PMAX=0.D0
        PMIN=1.D12
        DO NC=1,NC0
          VAL=VAL+PAS
          SR=0.D0
          SI=0.D0
          DO NT=1,npt
            SR=SR + spin(ksx,nt)*cos(nt*val) + spin(ksy,nt)*sin(nt*val) 
            SI=SI + spin(ksx,nt)*sin(nt*val) - spin(ksy,nt)*cos(nt*val) 
          ENDDO
          PP=SR*SR+SI*SI
          IF(PP.GT. PMAX) THEN
            PMAX=PP
            KMAX=NC
          ELSEIF(PP.LT. PMIN) THEN
            PMIN=PP
          ENDIF
          SPEC(NC)=PP
          Q(NC)=val/(2.d0*pi)
        ENDDO

        IF (PMAX .GT. PMIN) THEN
          IF (KMAX .LT. mxc) THEN
            DEC=0.5D0 * (SPEC(KMAX-1)-SPEC(KMAX+1))
     >      /(SPEC(KMAX-1) - 2.D0 *SPEC(KMAX)+SPEC(KMAX+1))
          ELSE
            DEC=0.5D0 
          ENDIF
          YNU= ANUI + (DBLE(KMAX)+DEC-0.5D0) * DELNU
        ELSE
           YNU = 0.D0
        ENDIF

      write(*,*) '---------------------------'
      write(*,*) '---------------------------'
      write(*,*) ' Number of channels : ',nc0
      write(*,*) ' Number of turns : ',npt
      write(*,*) ' Q/1-Q = ',ynu,'/',1.d0-ynu
      write(*,*) '---------------------------'
      write(*,*) '---------------------------'
      
      if(idluni(lunW)) 
     >   open(unit=lunW,file='spectrum.Out')
      write(lunW,fmt='(1p,2(e16.8,1x),i6)') (q(j),spec(j),j,j=1,nc0)
      close(lunW)

      Qs(jspec) = ynu
      write(lunW2,fmt='(i4,1p,2e14.6)') jspec, Qs(jspec)
      jspec = jspec+1
      ynu = 0.d0
      if(jspec.le.Nspec .and. .not.okend) goto 7

 22   continue
      write(*,*) ' Job completed !'
      close(lunW2)
      close(lunR)
      RETURN
      END
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
