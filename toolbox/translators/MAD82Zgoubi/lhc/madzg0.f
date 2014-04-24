C23456789012345678901234567890123456789012345678901234567890123456789012

      character fnr*80, fnw*80, txt*132
      parameter(lr=9, lw=10)
      parameter (mxl=15000)
      dimension iq(mxl)
      character kley*4, name*16
      parameter(mxk=15)
      character*4 kle(mxk)
      data kle/ 'DRIF', 'RBEN', 'SBEN', 'QUAD', 'SEXT', 'OCTU', 'MULT',
     >'SROT', 'YROT', 'MARK', 'KICK', 'HMON', 'VMON','HKIC', 'VKIC' /


C      write(*,*) ' name of the survey file to be translated:'
C      read(*,fmt='(A)') fnr
C      write(*,*) ' name of the zgoubi data file:'
C      read(*,fmt='(A)') fnw

      fnr='survey'
      fnw='zgoubi.dat'

      OPEN(UNIT=lr,FILE=fnr,status='old',err=97)
 80   OPEN(UNIT=lw,FILE=fnw,STATUS='NEW',IOSTAT=IOS)
      IF(IOS .NE. 0) THEN
        CLOSE(lw,STATUS='DELETE')
        OPEN(UNIT=lw,FILE=fnw,STATUS='OLD')
        CLOSE(lw,STATUS='DELETE')
        GOTO 80
      ENDIF
     
      read(lr,fmt='(A)',err=99) txt    
      read(lr,fmt='(A)',err=99) txt
      write(*,*) txt
      write(lw,fmt='(A60)') txt
      call objet(lw,bro)
      write(*,*) ' Fringe fields (y=1 / n=0) :'
      read(*,*) frf
      if( frf .ne. 1.) frf = 0.

      it = 1

      do i=1,4
        read(lr,fmt='(A)',err=99) txt
      enddo        

      write(*,*) ' Now translating. Wait...'

      ir = 0
      noel = 0
 86   continue
        read(lr,fmt='(A4,A16,F12.6,3E16.9,/,5E16.9)',end=85,err=99)
     >   kley, name, xl, ang, ak1, ak2, tilt, e1, e2, h1, h2
        read(lr,fmt='(1P,4E16.9)',end=85,err=99) bid,bid,bid,s   
!        write(*,fmt='(A4,A16,F12.6,3E16.9,/,5E16.9)')
!     >   kley, name, xl, ang, ak1, ak2, tilt, e1, e2, h1, h2
!        write(*,fmt='(E16.9)') s
        read(lr,fmt='(A)',end=85,err=99) txt
        ir = ir + 1

        DO 88 IK=1,MXK
          IF(KLEY .EQ. KLE(IK)) THEN
            NOEL = NOEL+1
            IF( NOEL .EQ. MXL+1) THEN
              TOMANY=.TRUE.
              NOEL=1
            ENDIF
            IQ(NOEL) =  IK
            GOTO 87
          ENDIF
 88     CONTINUE

        WRITE(*,100) KLEY
 100    FORMAT(/,10X,' PGM PRINCIPAL : ARRET SUR CLE   ',A,//,128(1H*))
        stop
 
 87     CONTINUE
        it = it + 1
        call lmnt
     >    (lw,bro,frf,ik,noel,kley,name,xl,ang,ak1,ak2,tilt,e1,e2,h1,h2,
     >                                                             it)

      goto 86

 85   continue
      write(lw,fmt='(A)') '''FAISCNL'''
      write(lw,fmt='(A)')   'zgoubi.fai'
      it = it + 1
      write(lw,fmt='(A)') '''MATRIX'''
      write(lw,fmt='(A)')   '1 11 '
      it = it + 1
      write(lw,fmt='(A)') '''END'''
      write(*,*) ' Read ',ir,' elements from the survey file '
      write(*,*) fnr
      write(*,*) ' Translated to',it,
     >' elements into the zgoubi data file'
      write(*,*) fnw
      stop ' end of survey file '
      
 97   stop ' error open survey file'
 99   stop ' error during read of survey file'
      end
        
      subroutine objet(lw,bro)
!      bro = 7. / 2.99792458 * 1.E4   ! T.m
      bro = sqrt(7000*(7000+2*.93827231))/.299792458    ! T.m -> 7 TeV proton
 1    continue
      write(*,*)  ' problem rigidity :'
      read(*,*,err=1) bro
      write(*,FMT='(F17.4,A)') bro,' T.m'
      write(lw,fmt='(A)') '''OBJET'''
      write(lw,FMT='(F17.4)') bro*1000.
      write(lw,fmt='(A)') '5' 
      write(lw,fmt='(A)') 
     >' .001 .01 .001 .01 0. .00001'
!     > '3.381847E-02 .005  1.714452E-02 .005 0. .00001'
      write(lw,fmt='(A)') '0. 0. 0. 0. 0. 1.'
      return
      end 

      subroutine lmnt
     >(lw,bro,frf,ik,noel,kley,name,xl,ang,ak1,ak2,tilt,e1,e2,h1,h2,
     >                                                         it)
      character name*16

      character*60 txfd, txfq
C----- Coeff champ de fuite Dipole LHC White book
!      data txfd /'6 -1.9067 6.0897 -2.7987 -.27635  -3.63676 7.3183' /
!      data txfq /'6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723'/
!      data txfd /'4  .1455   2.2670  -.6395  1.1558  0. 0.  0.'/ 
      data txfq 
     >/'6 -.010967  5.464823  .996848 1.568787 -5.671630 18.505734' /
      data txfd 
     >/'6 .015527 3.874961 -2.362230 2.978209 12.604429 15.025689' /
      parameter (pi=3.14159265389793238)
      parameter (cm=100., rd=180./pi, tkg=10.)
      parameter(i0=0, i1=1, i2=2, i4=4, i6=6)
      parameter(x0=0.,x1=1.,x2=2.,x3=3.,x4=4.,x5=5.,x6=6.,x7=7.,x8=8.)
      parameter(x10=10.,x20=20.,x999=999.)

      character txt*60

C----- 1    2    3    4    5    6    7    8    9    10   11   12   13
C----- DRIF RBEN SBEN QUAD SEXT OCTU MULT SROT YROT MARK KICK HMON VMON
C----- 14   15
C----- HKIC VKIC
!      goto (1,21,3,4,5,6,7,8,9,10,11,12,12,11,11) ik
       goto (1,2 ,3,4,5,6,7,8,9,10,11,12,12,11,11) ik

 1    continue
C----- DRIF
      if(xl .eq. 0) then
        it = it - 1
      else
        write(lw,fmt='(A,T12,A,T22,A)') '''DRIFT''',kley,name
        write(lw,fmt='(F12.4)') xl*cm
      endif
      goto 99
 2    continue
C----- RBEN -> MULTIPOL
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif
      write(lw,fmt='(A,T12,A,T22,A)') '''MULTIPOL''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .Dip'
      ro=xl/2./sin(ang/2.)
      b1 = bro / ro           ! dipole field (T)
!      b1 = ang * bro / xl 
      write(lw,fmt='(F8.2,F6.2,2F14.8,4F9.4)')
     >xl*cm,x10,b1*10., x0,x0,x0,x0,x0
      txt = '15.00 11.20  0.00 0.00 0.00 0.00 0.00' 
      if(frf .eq. 0.) txt = '0. 0. '//txt
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A)') '10.020E10  Dip'
      write(lw,fmt='(A)') '3 0. 0. 0.'
      goto 99    
 
 21   continue
C----- RBEN -> BEND
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif
      write(lw,fmt='(A,T12,A,T22,A)') '''BEND''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .Bend'
      ro=xl/2./sin(ang/2.)
      b=ang * bro / xl *10.  !bro/ro*tkg
      xxl=ro*ang*cm
      ro=ro*cm
      write(lw,fmt='(F14.7,F15.6,F15.8)') xxl,ro,b 
      te=ang/2.+e1
      write(lw,fmt='(2F6.2,F12.8)') x20,x8,te
      write(lw,fmt='(A)') txfd
      ts=ang/2.+e2
      write(lw,fmt='(2F6.2,F12.8)') x20,x8,ts
      write(lw,fmt='(A)') txfd
      step = ro * sin(ang/2.) * 2.
      write(lw,fmt='(F12.6,A4)') step / 10.,'  Bend'
      write(lw,fmt='(A)') '1 0. 0. 0.'
      goto 99

 3    continue
C----- SBEN
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif
      write(lw,fmt='(A,T12,A,T22,A)') '''BEND''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .Bend'
      ro=xl/ang
      b=bro/ro*tkg
      write(lw,fmt='(F14.7,2F10.4)') xl*cm,ro*cm,b 
      te=e1
      write(lw,fmt='(2F6.2,2F10.7)') x20,x8,te
      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
      ts=e2
      write(lw,fmt='(2F6.2,2F10.7)') x20,x8,ts
      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
      write(lw,fmt='(A)') '1.5  Bend'
      write(lw,fmt='(A)') '1 0. 0. 0.'
      goto 99
 4    continue
C----- QUAD
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif
      write(lw,fmt='(A,T12,A,T22,A)') '''MULTIPOL''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .Quad'
      write(lw,fmt='(F8.2,F6.2,2F14.8,4F9.4)')
     >xl*cm,x10,x0,ak1*bro,x0,x0,x0,x0
      txt = '8.0 5.60  1.00 0.00 0.00 0.00 0.00'
      if(frf .eq. 0.) txt = '0. 0. '//txt
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfq
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfq
      write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A)') '10.020E10  Quad'
      write(lw,fmt='(A)') '1 0. 0. 0.'
      goto 99
 5    continue
C----- SEXT ***
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif
      write(lw,fmt='(A,T12,A,T22,A)') '''MULTIPOL''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .Sext'
      write(lw,fmt='(F8.2,F6.2,F14.8,F12.6,4F9.4)')
     >xl*cm,x10,x0,x0,x0,x0,x0,x0
      write(lw,fmt='(A)') '0.00 0.00  1.00 0.00 0.00 0.00 0.00'
      write(lw,fmt='(A)') txfq
      write(lw,fmt='(A)') '0.00 0.00  1.00 0.00 0.00 0.00 0.00'
      write(lw,fmt='(A)') txfq
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A)') '10.020E10  Sext'
      write(lw,fmt='(A)') '1 0. 0. 0.'
      goto 99
 6    continue
C----- OCTU ***
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif
      write(lw,fmt='(A,T12,A,T22,A)') '''MULTIPOL''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .Octu'
      write(lw,fmt='(F8.2,F6.2,F14.8,F12.6,4F9.4)')
     >xl*cm,x10,x0,x0,x0,x0,x0,x0
      write(lw,fmt='(A)') '0.000 0.000  1.00 0.00 0.00 0.00 0.00'
      write(lw,fmt='(A)') txfq
      write(lw,fmt='(A)') '0.000 0.000  1.00 0.00 0.00 0.00 0.00'
      write(lw,fmt='(A)') txfq
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A)') '10.020E10  Octu'
      write(lw,fmt='(A)') '1 0. 0. 0.'
       goto 99
 7    continue
C----- MULT
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif
      write(lw,fmt='(A,T12,A,T22,A)') '''MULTIPOL''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .Mult'
      write(lw,fmt='(F8.2,F6.2,2F14.8,4F9.4)')
     >xl*cm,x10,x0,ak1*bro,x0,x0,x0,x0
      write(lw,fmt='(A)') '00.00 00.00  1.00 0.00 0.00 0.00 0.00'
      write(lw,fmt='(A)') txfq
      write(lw,fmt='(A)') '00.00 00.00  1.00 0.00 0.00 0.00 0.00'
      write(lw,fmt='(A)') txfq
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A)') '10.020E10  Mult'
      write(lw,fmt='(A)') '1 0. 0. 0.'
       goto 99
 8    continue
C----- SROT ***
      write(lw,fmt='(A,T12,A,T22,A)') '''TRAROT''',kley,name
      goto 99
 9    continue
C----- YROT ***
      write(lw,fmt='(A,T12,A,T22,A)') '''TRAROT''',kley,name
      goto 99
 10   continue
C----- MARK ***
      it = it - 1
      goto 99
 11   continue
C----- KICK ***
      xxl = xl
      if(xl .eq. 0) xxl = 1.e-6
      if(ik.eq.14) then
c--------- hkicker
        diptlt = 0.
        b = -bro*tilt/xxl*tkg
      elseif(ik.eq.15) then
c--------- vkicker
        diptlt = pi/2.
        b = -bro*e1/xxl*tkg
      endif
      write(lw,fmt='(A,T12,A,T22,A)') '''MULTIPOL''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .kicker'
      write(lw,fmt='(E12.4,F6.2,E14.6,F12.6,4F9.4)')
     >xxl*cm,x10,b,x0,x0,x0,x0,x0
      txt = '.0 .0  1.00 0.00 0.00 0.00 0.00'
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(f12.9,A)') diptlt,' 0. 0. 0. 0. 0.'
      write(lw,fmt='(A)') '10.020E10  Kick'
      write(lw,fmt='(A)') '1 0. 0. 0.'
      goto 99
 12   continue
C----- HMON, VMON
      write(lw,fmt='(A,T12,A,T22,A)') '''DRIFT''',kley,name
      write(lw,fmt='(F12.4)') xl*cm
      goto 99

 99   return
      end
