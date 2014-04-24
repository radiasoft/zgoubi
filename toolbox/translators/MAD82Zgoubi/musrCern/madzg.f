C23456789012345678901234567890123456789012345678901234567890123456789012
      implicit double precision (a-h,o-z)
      character fnr*80, fnw*80, txt*132
      parameter(lr=9, lw=10)
      parameter (mxl=15000)
      dimension iq(mxl)
      character kley*4, name*16
      parameter(mxk=15)
      character kle(mxk)*4, ny*1
      logical TOMANY

      data kle/ 'DRIF', 'RBEN', 'SBEN', 'QUAD', 'SEXT', 'OCTU', 'MULT',
     >'SROT', 'YROT', 'MARK', 'KICK', 'HMON', 'VMON','HKIC', 'VKIC' /

C      write(*,*) ' name of the survey file to be translated:'
C      read(*,fmt='(A)') fnr
C      write(*,*) ' name of the zgoubi data file:'
C      read(*,fmt='(A)') fnw

      fnr='survey'
      fnw='trad.out'

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

      write(*,*) 
      write(*,*) ' Fringe field coeff :',
     >             ' lhc/recycler/musr/default (1/2/3/99)'
 2    read(*,*,err=2) kff
C      if( kff .ne. 1 .and. kff .ne. 2) kff = 3
      write(*,*) ' Fringe field coeff : option set to ',kff
      call lmntff(kff)

      write(*,*) 
      write(*,*) ' Fringe fields on/off/default (1/0/99) :'
 3    read(*,*,err=3) frf
      if( frf .ne. 1.d0) frf = 0.d0
      if(frf.eq.1.d0) write(*,*) ' Fringe fields on '
      if(frf.eq.0.d0) write(*,*) ' Fringe fields off '

      write(*,*) 
      write(*,*) ' Option KPOS 1/2/3/default (1/2/3/99) :'
 4    read(*,*,err=4) kpos
      if( kpos .lt. 1 .or. kpos .gt. 3) kpos=2
      write(*,*) ' Option KPOS set to ',kpos
      xce=0.d0
      yce=0.d0
      ale=0.d0
      if( kpos.eq.2) then
 5      write(*,*) ' Give 3 alignement values xce,yce,ale :'
        read(*,*,err=5) xce,yce,ale
      endif
      call lmnkpo(xce,yce,ale,kpos)

      write(*,*)
      write(*,*) ' Convert DRIFT to zero-field MULTIPOLE (n/y) :'
      read(*,*,err=6) ny
 6    continue
      if(ny.eq.'Y') ny='y' 
      if(ny.eq.'N') ny='n' 
      call lmndri(ny)

      it = 1

      do i=1,4
        read(lr,fmt='(A)',err=99) txt
      enddo        

      write(*,*) 
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
C 100    FORMAT(/,10X,' PGM PRINCIPAL : ARRET SUR CLE   ',A,//,128(1H*))
 100    FORMAT(/,10X,' Key ',A,' not translated...')
        goto 86 
C        stop
 
 87     CONTINUE
        it = it + 1
        call lmnt
     >    (lw,bro,frf,ik,noel,kley,name,xl,ang,ak1,ak2,tilt,e1,e2,h1,h2,
     >                                      it)

      goto 86

 85   continue
      write(lw,fmt='(A)') '''FAISCNL'''
      write(lw,fmt='(A)')   'zgoubi.fai'
      it = it + 1
      write(lw,fmt='(A)') '''MATRIX'''
      write(lw,fmt='(A)')   '2  11 '
      it = it + 1
      write(lw,fmt='(A)') '''END'''
      it = it + 1
      write(lw,fmt='(///,A)') '''REBELOTE'''
      write(lw,fmt='(A)')   '600 0.1 99'
      write(*,*) ' Read ',ir,' elements from the survey file ',fnr
      write(*,*) ' Translated to',it,
     >' elements into the zgoubi data file "trad.out"', fnw
      stop ' end of survey file '
      
 97   stop ' error open survey file'
 99   stop ' error during read of survey file'
      end
        
      subroutine objet(lw,bro)
      implicit double precision (a-h,o-z)
     
 1    continue

c      write(*,*) '  Give kinetic energy of protons  (GeV) '
c      read(*,*,err=1) T
c      bro = sqrt(T*(T+2.d0*.93827231d0))/.299792458d0    ! T.m 
C      write(*,*) '  Give Brho (T.m) '
C      read(*,*,err=1) bro 

!        bro = 29.65014391d0
        bro = 1504.164917d0  ! LHC, 450 GeV

      write(*,*)  ' problem rigidity :'
      write(*,FMT='(F19.6,A,/,A)') 
     >      bro,' T.m','  Enter desired value'
      read(*,*,err=1) bro
   
      write(lw,fmt='(A)') '''OBJET'''
      write(lw,FMT='(F19.6)') bro*1000.d0
      write(lw,fmt='(A)') '5' 
      write(lw,fmt='(A)') '.001 .001 .001 .001 0. .000001  '
      write(lw,fmt='(A)') '0. 0. 0. 0. 0. 1.'
      return
      end 

      subroutine lmnt
     >(lw,bro,frf,ik,noel,kley,name,xl,ang,ak1,ak2,tilt,e1,e2,h1,h2,
     >                                                         it)
      implicit double precision (a-h,o-z)
      character name*16, ny*1

      character*80 txfd, txfq,fqlhc , fdlhc ,fqrec , fdrec,fd,fq
      character*80 txffd,txffq,ffqlhc,ffdlhc,ffqrec,ffdrec,ffd,ffq
      character*80 ffdmu, fdmu, ffqmua,ffqmus,fqmua
      save txfd,txfq,txffd, txffq
      logical drimul

      data ffq / '6.00  3.00 1.00 0.00 0.00 0.00 0.00 0. 0. 0. 0. ' / 
      data fq / '6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723'/
      data ffd / '10.00  4.0  0.800 0.00 0.00 0.00 0.00 0. 0. 0. 0. ' / 
      data fd /'4  .1455   2.2670  -.6395  1.1558  0. 0.  0.'/ 

!----------- Recycler fringe fields -----------------------------------------
!  quadrupole
      data ffdrec / '8.00  5.0  1.000 1.00 0.00 0.00 0.00 0. 0. 0. 0. '/ 
      data fdrec /'4  0.09650  3.76444 -0.70378  1.31734  0. 0. 0.' / 
!  dipole
      data ffqrec / '6.00  3.00 1.00 0.00 0.00 0.00 0.00 0. 0. 0. 0. ' / 
      data fqrec / '6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723' / 

!----------- LHC fringe fields -----------------------------------------
!  dipole White book
      data ffdlhc /'15.00 11.20  0.00 0.00 0.00 0.00 0.00 0. 0. 0. 0.' / 
      data fdlhc 
     >/'6 .015527 3.874961 -2.362230 2.978209 12.604429 15.025689' /
!  quadrupole Saclay
      data ffqlhc / '8.0 5.60  1.00 0.00 0.00 0.00 0.00 0. 0. 0. 0. ' / 
      data fqlhc 
     >/'6 -.010967  5.464823  .996848 1.568787 -5.671630 18.505734' /
!----------- LHC fringe fields -----------------------------------------

!----------- muon collider fringe fields -------------------------------
!  dipoles
      data ffdmu /'15.00 11.20  0.00 0.00 0.00 0.00 0.00 0. 0. 0. 0.' / 
      data fdmu 
     >/'6 .015527 3.874961 -2.362230 2.978209 12.604429 15.025689' /
!  arc quadrupoles (length=0.5, radius=31mm)
      data ffqmua / '9.0 9.0  1.00 0.00 0.00 0.00 0.00 0. 0. 0. 0.' / 
      data fqmua 
     >/'6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723  arc' /
!  LSS and matching cell quadrupoles (length=0.35, radius=89mm)
      data ffqmus / '17.0 20.0  1.00 0.00 0.00 0.00 0.00 0. 0. 0. 0.' / 
!----------- muon collider fringe fields -------------------------------

      parameter (cm=100.d0, tkg=10.d0)
      parameter(i0=0, i1=1, i2=2, i4=4, i6=6)
      parameter(x0=0.d0,x1=1.d0,x2=2.d0,x3=3.d0,x4=4.d0,x5=5.d0,x6=6.d0)
      parameter(x7=7.d0,x8=8.d0,x10=10.d0,x20=20.d0,x999=999.d0)

      character txt*80

      pi = 4.d0*atan(1.d0)
      zero=0.d0

C----- 1    2    3    4    5    6    7    8    9    10   11   12   13
C----- DRIF RBEN SBEN QUAD SEXT OCTU MULT SROT YROT MARK KICK HMON VMON
C----- 14   15
C----- HKIC VKIC
!      goto (1,21,3,4,5,6,7,8,9,10,11,12,12,11,11) ik
       goto (1,2 ,3,4,5,6,7,8,9,10,11,12,12,11,11) ik

 1    continue
C----- DRIF      
      if(xl .eq. 0.d0) then
        it = it - 1
      else
        if(drimul) then
          dum=1.d-20
          write(lw,fmt='(A,T12,A,T22,A)') '''MULTIPOL''',kley,name
          write(lw,fmt='(I1,A)') i0,'  .Drift'
          write(lw,fmt='(G10.4,G8.2,G16.10,9G7.1)')
     >    xl*cm,x10,x0,dum,x0,x0,x0,x0,x0,x0,x0,x0
          txt = '0. 0. '//txffq
          write(lw,fmt='(A)') txt
          write(lw,fmt='(A)') txfq
          write(lw,fmt='(A)') txt
          write(lw,fmt='(A)') txfq
          write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
          write(lw,fmt='(A)') '3.003E10  Drift'
          write(lw,fmt='(A)') '1 0. 0. 0.'
        else
          write(lw,fmt='(A,T12,A,T22,A)') '''DRIFT''',kley,name
          write(lw,fmt='(F12.4)') xl*cm
        endif
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

!      xlmag = xl
!      ro=xlmag/(2.d0*sin(ang/2.d0))
!      xlarc = ro * ang
      xlarc = xl
      ro = xlarc /ang
      xlmag = 2.d0 * ro * sin (ang/2.d0)

      b1 = bro / ro                              ! dipole field (T)
! Takes care of the offset due to grad in combined fcntn magnets
      if( ak1 .ne. 0.d0) then
        sqrk = sqrt(abs(ak1))
        sqrkl = sqrk * xlarc
        if(ak1 .ge. 0.d0 ) then
          si = sin(sqrkl)
          co = cos(sqrkl)
        else
          si = sinh(sqrkl)
          co = cosh(sqrkl)
        endif
        off = ( si/2.d0/sqrk/(1.d0-co)- 1.d0/xlarc/ak1)*ang
!        b1 = b1 + ak1 * off * Bro
!        b1 =  ak1 *Bro * si/2.d0/sqrk/(1.d0-co)*ang
!        write(*,fmt='(a,4f16.10)') name,off,ro,bro / ro, b1
      endif
! quadrupole field at .1 m (T)
      b2 = ak1*bro * x10/100.                     
! sextupole field at .1 m (T)
      b3 = ak2 * bro *(x10/100.)**2 / 2.d0
C         write(*,*) ak2, bro, b3, x10

      write(lw,fmt='(F9.4,F7.2,3F13.8,7F4.1)')
     >xlmag*cm, x10, b1*10.d0, b2*10.d0, b3*10.d0, x0,x0,x0,x0,x0,x0,x0

      txt=txffd
      if(frf .eq. 0.) txt='0. 0. '//txt
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A,2X,A)') ' 10.010E10    Dip',name
      if(kpos.eq.3) then
        write(lw,fmt='(I,A,1P,2G18.10)') kpos,' 0. ',-off*100.,zero
C        write(lw,fmt='(I,A,1P,2G18.10)') kpos,' 0. ',-off*100.,ang/2.
      else
        write(lw,fmt='(I,1P,3G18.10)') kpos,xce,yce,ale
      endif
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
C      xxl=ro*ang*cm
      ro=ro*cm

      write(lw,fmt='(2F14.7,F15.8)') xl*cm,tilt,b 
      te=ang/2.+e1
      write(lw,fmt='(2F6.2,F12.8)') x20,x8,te
      write(lw,fmt='(A)') txfd
      ts=ang/2.+e2
      write(lw,fmt='(2F6.2,F12.8)') x20,x8,ts
      write(lw,fmt='(A)') txfd
C      step = ro * sin(ang/2.) * 2.
C      write(lw,fmt='(F12.6,A4)') step / 10.,'  Bend'
      write(lw,fmt='(A,2X,A)') ' 10.010E10   Bend',name
      if(kpos.eq.3) then
        write(lw,fmt='(A)') '3 0. 0. 0.'
      else
        write(lw,fmt='(I,1P,3G18.10)') kpos,xce,yce,ale
      endif
      goto 99

 3    continue
C----- SBEN -> BEND
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif

      write(lw,fmt='(A,T12,A,T22,A)') '''BEND''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .Bend'
      ro=xl/ang
      xxl = 2.d0*ro*sin(ang/2.d0)
      b=bro/ro*tkg
              write(*,*) ' ********* tilt ',tilt
      write(lw,fmt='(2F14.7,F10.4)') xxl*cm,tilt,b 
      te=e1
      write(lw,fmt='(2F6.2,2F10.7)') x20,x8,te
      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
      ts=e2
      write(lw,fmt='(2F6.2,2F10.7)') x20,x8,ts
      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
      write(lw,fmt='(A,2X,A)') '1.5  Bend',name
      if(kpos.eq.3) then
        write(lw,fmt='(A)') '3 0. 0. 0.'
      else
        write(lw,fmt='(I,1P,3G18.10)') kpos,xce,yce,ale
      endif
      goto 99
 4    continue
C----- QUAD
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif
      write(lw,fmt='(A,T12,A,T22,A)') '''MULTIPOL''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .Quad'
      write(lw,fmt='(F10.4,F6.2,2F16.10,8F4.1)')
     >xl*cm,x10,x0,ak1*bro,x0,x0,x0,x0,x0,x0,x0,x0
      txt = txffq
      txt=txt//name
      if(frf .eq. 0.) txt = '0. 0. '//txt
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfq
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfq
      write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A,2X,A)') '10.010E10  Quad',name
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
      b3 = ak2 * bro *(x10/100.)**2 / 2.d0
      write(*,*) ' SEXT         b3 = ',b3,' kG'
      write(lw,fmt='(F10.4,F6.2,2F10.4,F16.10,7F4.1)')
     >xl*cm,x10,x0,x0,b3,x0,x0,x0,x0,x0,x0,x0
      write(lw,fmt='(A)') '0.00 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') '0.00 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A,2X,A)') '10.010E10  Sext',name
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
      write(lw,fmt='(F10.4,F6.2,2F16.10,4F9.4)')
     >xl*cm,x10,x0,x0,x0,x0,x0,x0
      write(lw,fmt='(A)') '0.000 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') '0.000 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A,2X,A)') '20.020E10  Octu',name
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
      write(lw,fmt='(F10.4,F6.2,2F16.10,8F4.1)')
     >xl*cm,x10,x0,ak1*bro,x0,x0,x0,x0
      write(lw,fmt='(A)') '00.00 00.0  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') '00.00 00.0  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
C      write(lw,fmt='(A,2X,A)') '20.020E10  Mult',name
      ten=10.d0
      write(lw,fmt='(F8.2)') ten
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
        diptlt = pi/2.d0
        b = -bro*e1/xxl*tkg
      endif
      write(lw,fmt='(A,T12,A,T22,A)') '''MULTIPOL''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .kicker'
      write(lw,fmt='(E12.4,F6.2,E14.6,F12.6,8F4.1)')
     >xxl*cm,x10,b,x0,x0,x0,x0,x0,x0,x0,x0,x0
      txt = '.0 .0  1.00 0.00 0.00 0.00 0.00 0. 0. 0. 0.'
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(f12.9,A)') diptlt,' 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A)') '20.020E10  Kick'
      write(lw,fmt='(A)') '1 0. 0. 0.'
      goto 99
 12   continue
C----- HMON, VMON
      write(lw,fmt='(A,T12,A,T22,A)') '''DRIFT''',kley,name
      write(lw,fmt='(F12.4)') xl*cm
      goto 99
 99             write(6,fmt='(A,T12,A,)') kley,name
        return

      entry lmntff(kff)

      if(kff .eq.1) then
!-------- LHC
        txffd = ffdlhc
        txfd = fdlhc
        txffq = ffqlhc
        txfq = fqlhc        
      else if(kff .eq.2) then
!-------- Recycler      
        txffd = ffdrec
        txfd = fdrec
        txffq = ffqrec
        txfq = fqrec
      else if(kff .eq.3) then
!-------- muon collider 
        txffq = ffqmua
        txfq = fqmua
        txffd = ffdlhc
        txfd = fdlhc
      else
        txffd = ffd
        txfd = fd
        txffq = ffq
        txfq = fq
      endif

      write(*,fmt='(A)') txffd ,txfd,txffq ,txfq

      return

      entry lmnkpo(xcei,ycei,alei,kposi)
      kpos=kposi
      xce=xcei
      yce=ycei
      ale=alei
      return

      entry lmndri(ny)
      drimul=ny.eq.'y'
      return

      end
