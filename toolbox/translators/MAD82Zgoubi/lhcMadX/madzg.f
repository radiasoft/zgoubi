C23456789012345678901234567890123456789012345678901234567890123456789012
      implicit double precision (a-h,o-z)
      character fnr*80, fnw*80, txt*132
      parameter(lr=9, lw=10)
      parameter (mxl=15000)
      dimension iq(mxl)
      character kley*4, name*16
      character*17 txa17, txb17
      parameter(mxk=19)
      character kle(mxk)*4, ny*1
      logical TOMANY

      character txt262*262

      character nydri,nysex,nykic,nycon 
      save kff,frf,kpos,nydri,nysex, kbm,nykic,nycon 
      data kff,frf,kpos,nydri,nysex,kbm,nykic,nycon 
     >   / 1, 0.d0, 3, 'n', 'n', 2, 'y', 'y' /

      data xce,yce,ale /0.d0, 0.d0, 0.d0/

      data kle/ 'DRIF', 'RBEN', 'SBEN', 'QUAD', 'SEXT', 'OCTU', 'MULT',
     >'SROT', 'YROT', 'MARK', 'KICK', 'HMON', 'VMON','HKIC', 'VKIC',
     >'MONI', 'RCOL', 'ECOL' , 'RFCA'/

C      write(*,*) ' name of the survey file to be translated:'
C      read(*,fmt='(A)') fnr
C      write(*,*) ' name of the zgoubi data file:'
C      read(*,fmt='(A)') fnw

      write(*,*) 
 22   write(*,*) ' Input data file name : survey/madxData (1/2) : ' 
      read(*,*,err=22) irep
      if(irep.eq.1) then 
        fnr='survey'
      else
        fnr='madxData'
      endif
      fnw='trad.out'

      OPEN(UNIT=lr,FILE=fnr,status='old',err=97)
 850  OPEN(UNIT=lw,FILE=fnw,STATUS='NEW',IOSTAT=IOS)
      IF(IOS .NE. 0) THEN
        CLOSE(lw,STATUS='DELETE')
        OPEN(UNIT=lw,FILE=fnw,STATUS='OLD')
        CLOSE(lw,STATUS='DELETE')
        GOTO 850
      ENDIF
     
      if(fnr.eq.'survey') then
        read(lr,fmt='(A)',err=99) txt    
        read(lr,fmt='(A)',err=99) txt
        write(*,*) txt
      elseif(fnr.eq.'madxData') then
        read(lr,fmt='(A)',err=99) txt    
        write(*,*) txt(1:80)
        read(lr,fmt='(A)',err=99) txt
        write(*,*) txt(1:80)
      endif

C----- Write title = first line in zgoubi data file
      write(lw,fmt='(A60)') txt
      call objet(lw,bro)

      call lmntff(kff)
      call lmnkpo(xce,yce,ale,kpos)
      call lmndri(nydri)
      fb3 = 1.d0
      if(nysex.eq.'y') fb3=0.d0
      call lmnfb3(fb3)
      call lmnkbm(kbm)
      cki = 1.d0
      if(nykic.eq.'y') cki=0.d0
      call lmnkic(cki)
      call lmncon(nycon)

C----- Options
 888  continue
      WRITE(6,100) kff,frf,kpos,nydri,nysex,kbm,nykic,nycon
 100  FORMAT(//,3X,75(1H*),//,20X,' MENU - Graphic processing :' ,//
     > ,9X,'        ' ,/
     1 ,9X,' 2   Fringe field ',
     >            '- lhc/recycler/musr/default (1/2/3/99) : ',T70,I1,/
     2 ,9X,' 3   Fringe fields on/off/default (1/0/99) : ',T70,F2.0 ,/
     3 ,9X,' 4   Option KPOS 1/2/3/default (1/2/3/99) : ',T70,I1,/
     4 ,9X,' 5   Convert DRIFT to zero-field MULTIPOLE (n/y) : ',T70,A,/
     5 ,9X,' 6   All sextupoles set to zero (n/y) : ',T70,A,/
     7 ,9X,' 7   Translate SBEND to BEND or to MULTIPOLE (1/2) : ',
     >                                                       T70,I1,/
     8 ,9X,' 8   All kickers set to zero (n/y) : ',T70,A,/
     9 ,9X,' 9   Concatenate drifts (y/n) :  ',T70,A,/
     2 ,3X,75(1H*),/)

      WRITE(6,101) ' * Option  number ("return" if ok) : '
 101  FORMAT(A,$)
      READ(5,201,ERR=888) IOPT
 201  FORMAT(I2)
      GOTO (777,20, 30, 40, 50, 60, 70, 80, 90) IOPT
      GOTO 777


 20   write(*,*) 
      write(*,*) ' Fringe field coeff :',
     >             ' lhc/recycler/musr/default (1/2/3/99)'
 2    read(*,*,err=2) kff
      if( kff .ne. 1 .and. kff .ne. 2) kff = 3
      write(*,*) ' Fringe field coeff : option set to ',kff
      call lmntff(kff)
      GOTO 888

 30   write(*,*) 
      write(*,*) ' Fringe fields on/off/default (1/0/99) :'
 3    read(*,*,err=3) frf
      if( frf .ne. 1.d0) frf = 0.d0
      if(frf.eq.1.d0) write(*,*) ' Fringe fields on '
      if(frf.eq.0.d0) write(*,*) ' Fringe fields off '
      GOTO 888

 40   write(*,*) 
      write(*,*) ' Option KPOS 1/2/3/default (1/2/3/99) :'
 4    read(*,*,err=4) kpos
      if( kpos .lt. 1 .or. kpos .gt. 3) kpos=2
      write(*,*) ' Option KPOS set to ',kpos
      xce=0.d0
      yce=0.d0
      ale=0.d0
      if( kpos.eq.2) then
 41      write(*,*) ' Give 3 alignement values xce,yce,ale :'
        read(*,*,err=41) xce,yce,ale
      endif
      call lmnkpo(xce,yce,ale,kpos)
      GOTO 888

 50   write(*,*)
      write(*,*) ' Convert DRIFT to zero-field MULTIPOLE (n/y) :'
      read(*,*,err=5) nydri
 5    continue
      if(nydri.eq.'Y') nydri='y' 
      if(nydri.ne.'y') nydri='n' 
      write(*,*) '               option set to ',nydri
      call lmndri(nydri)
      GOTO 888

 60   write(*,*)
      write(*,*) ' All sextupoles set to zero (n/y) :'
      read(*,*,err=6) nysex
 6    continue
      if(nysex.eq.'Y') nysex='y' 
      if(nysex.ne.'y') nysex='n' 
      write(*,*) '               option set to ',nysex
      fb3 = 1.d0
      if(nysex.eq.'y') fb3=0.d0
      call lmnfb3(fb3)
      GOTO 888

 70   write(*,*)
      write(*,*) ' Translate SBEND to BEND or to MULTIPOLE (1/2) :'
      read(*,*,err=7) kbm
 7    continue
      if( kbm .ne. 2   ) kbm = 1
      if(kbm.eq.1   ) write(*,*) ' BEND will be used for SBEND'
      if(kbm.eq.2) write(*,*) ' MULTIPOLE will be used for SBEND'
      call lmnkbm(kbm)
      GOTO 888

 80   write(*,*)
      write(*,*) ' All kickers set to zero (n/y) :'
      read(*,*,err=8) nykic
 8    continue
      if(nykic.eq.'Y') nykic='y' 
      if(nykic.ne.'y') nykic='n' 
      write(*,*) '               option set to ',nykic
      cki = 1.d0
      if(nykic.eq.'y') cki=0.d0
      call lmnkic(cki)
      GOTO 888

 90   write(*,*)
      write(*,*) ' Concatenate drifts (y/n) :'
      read(*,*,err=9) nycon
 9    continue
      if(nycon.eq.'N') nycon='n' 
      if(ny.ne.'n') ny='y' 
      write(*,*) '               option set to ',nycon
      call lmncon(nycon)
      GOTO 888

C---------------------------------------------
 777  continue
      it = 0

      do i=1,4
        read(lr,fmt='(A)',err=99) txt
      enddo        

      write(*,*) 
      write(*,*) ' Now translating. Wait...'

      ir = 0
      noel = 0
 86   continue
      if(fnr.eq.'survey') then 
        read(lr,fmt='(A4,A16,F12.6,3E16.9,/,5E16.9)',end=85,err=99)
     >  kley, name, xl, ang, ak1, ak2, tilt, e1, e2, h1, h2
        read(lr,fmt='(1P,4E16.9)',end=85,err=99) bid,bid,bid,s   
        read(lr,fmt='(A)',end=85,err=99) txt
!        write(*,fmt='(A4,A16,F12.6,3E16.9,/,5E16.9)')
!     >   kley, name, xl, ang, ak1, ak2, tilt, e1, e2, h1, h2
!        write(*,fmt='(E16.9)') s
      elseif(fnr.eq.'madxData') then
        read(lr,fmt='(A262)',end=85,err=85) txt262
c           write(*,*) txt262
        read(txt262,fmt='(2A17)',end=85,err=85) txa17,txb17
        kley = txa17(1:4)
        name = txb17(1:4)
        txt262=txt262(35:262)
        read(txt262,*,end=85,err=85)
     >  xl,ang,akl0,akl1,akl2,akl3,tilt,e1,e2,h1,h2,s
c        write(*,*) xl,ang,akl0,akl1,akl2,akl3,tilt,e1,e2,h1,h2,s
c        pause
c        read(lr,fmt='(2A17,12E19.0)',end=85,err=85)
c     >  txa17,txb17, xl,ang,akl0,akl1,akl2,akl3,tilt,e1,e2,h1,h2,s
c        kley = txa17(1:4)
c        name = txb17(1:4)
c        if(name.eq.'TCL.') then
c            write(*,*) '.......................',xl
c              stop
c         endif
c        if(name.eq.'MB.B') then
c              stop
c         endif
        if(xl.ne.0.d0) then
          ak0=akl0 / xl
          ak1=akl1 / xl
          ak2=akl2 / xl
          ak3=akl3 / xl
        else
        endif
      endif
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

        WRITE(*,102) KLEY
 102    FORMAT(/,10X,' Key ',A,' not translated...')
        goto 86 
C        stop
 
 87     CONTINUE
        it = it + 1
        call lmnt
     >    (lw,bro,frf,ik,noel,kley,name,xl,ang,ak1,ak2,tilt,e1,e2,h1,h2,
     >                                      ss,it)
        write(88,fmt='(A15,I6,2X,2A4,I6,3F10.2)')
     >                '  read lmnt # ',ir, kley,name,it,xl,ss,s

      goto 86

 85   continue
      write(lw,fmt='(A)') '''MATRIX'''
      write(lw,fmt='(A)')   '2  11 '
      it = it + 1

      write(lw,fmt='(A)') '''END'''
      it = it + 1

      write(lw,fmt='(///,A)') '''REBELOTE'''
      write(lw,fmt='(A)')   '399 0.1 99'
      write(lw,fmt='(A)')   'END'
      write(*,*) ' Read ',ir,' elements from the data file ',fnr
      write(*,*) ' Translated into',it,
     >' elements into the zgoubi data file ', fnw
      stop '  ///  Translation ended  ///         '
      
 97   stop ' error open input data file'
 99   stop ' error during read of input data file'
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

      broo = bro
      write(*,*)  ' problem rigidity :'
      write(*,FMT='(F19.6,A,/,A)') 
     >      bro,' T.m','  Enter desired value'
      read(*,*,end=2,err=2) bro
      goto 3
 2    continue
      bro=broo

 3    continue
      write(lw,fmt='(A)') '''OBJET'''
      write(lw,FMT='(F19.6)') bro*1000.d0
      write(lw,fmt='(A)') '5' 
      write(lw,fmt='(A)') '.0001 .001 .0001 .001 0. .000001  '
      write(lw,fmt='(A)') '0. 0. 0. 0. 0. 1.'

      write(lw,fmt='(A)') '''FAISCNL'''
      write(lw,fmt='(A)')   'b_zgoubi.fai'

      return
      end 

      subroutine lmnt
     >(lw,bro,frf,ik,noel,kley,name,xl,ang,ak1,ak2,tilt,e1,e2,h1,h2,
     >                                                        ss,it)
      implicit double precision (a-h,o-z)
      character kley*(*), name*16, ny*1

      character*80 txfd, txfq,fqlhc , fdlhc ,fqrec , fdrec,fd,fq
      character*80 txffd,txffq,ffqlhc,ffdlhc,ffqrec,ffdrec,ffd,ffq
      character*80 ffdmu, fdmu, ffqmua,ffqmus,fqmua
      character*20 txtsD,txtsQ, stpD, stpQ, sDlhc, sQlhc

      save txfd,txfq,txffd, txffq, stpD, stpQ, sDlhc, sQlhc
      logical drimul
      save drimul

      data ffq / '6.00  3.00 1.00 0.00 0.00 0.00 0.00 0. 0. 0. 0. ' / 
      data fq / '6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723'/
      data ffd / '10.00  4.0  0.800 0.00 0.00 0.00 0.00 0. 0. 0. 0. ' / 
      data fd /'4  .1455   2.2670  -.6395  1.1558  0. 0.  0.'/ 
      data stpD, stpQ / ' #20|60|20    Dip  ', ' #20|60|20    Quad ' /
      data txtsD, txtsQ / ' #20|60|20    Dip  ', ' #20|60|20    Quad ' /

!----------- Recycler fringe fields -----------------------------------------
!  quadrupole
      data ffdrec / '8.00  5.0  1.000 1.00 0.00 0.00 0.00 0. 0. 0. 0. '/ 
      data fdrec /'4  0.09650  3.76444 -0.70378  1.31734  0. 0. 0.' / 
!  dipole
      data ffqrec / '6.00  3.00 1.00 0.00 0.00 0.00 0.00 0. 0. 0. 0. ' / 
      data fqrec / '6  .1122 6.2671 -1.4982 3.5882 -2.1209 1.723' / 

!----------- LHC fringe fields -----------------------------------------
!  dipole White book
      data ffdlhc /'17.00 11.20  0.00 0.00 0.00 0.00 0.00 0. 0. 0. 0.' / 
      data fdlhc 
     >/'6 .015527 3.874961 -2.362230 2.978209 12.604429 15.025689' /
!  quadrupole Saclay
      data ffqlhc / '10.0 5.60  1.00 0.00 0.00 0.00 0.00 0. 0. 0. 0. ' / 
      data fqlhc 
     >/'6 -.010967  5.464823  .996848 1.568787 -5.671630 18.505734' /
      data sDlhc, sQlhc / '#30|160|30    Dip  ', '#30|160|30    Quad  '/
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
      save kpos, kbm, cik, fb3

      logical condri
      data condri /.false./
      save condri

      character*4 prec
      data prec, xlprec  / 'none' , 0.d0/ 
      save prec

      pi = 4.d0*atan(1.d0)
      zero=0.d0

C----- 1    2    3    4    5    6    7    8    9    10   11   12   13
C----- DRIF RBEN SBEN QUAD SEXT OCTU MULT SROT YROT MARK KICK HMON VMON
C----- 14   15   16   17   18   19
C----- HKIC VKIC MONI RCOL ECOL RFCA
       goto (1,2,2,4,5,6,7,8,9,10,11,12,12,11,11,12,13,13,1 ) ik

 1    continue
C----- DRIF      
      if(xl .eq. 0.d0) then
        it = it - 1
      else
        ss = ss+xl
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
c          write(lw,fmt='(A)') '3.003E10  Drift'
          write(lw,fmt='(A)') '#3|3|3  Drift'
          write(lw,fmt='(A)') '1 0. 0. 0.'
        else
          if(prec.eq. 'DRIF' .and. condri) then
            it = it-1
            xlprec = xlprec + xl
            backspace(lw)
          else
            write(lw,fmt='(A,T12,A,T22,A)') '''DRIFT''',kley,name
            prec = 'DRIF'
            xlprec = xl
          endif
          write(lw,fmt='(F14.6)') xlprec*cm
        endif
      endif
      goto 99

 2    continue
      if(kbm.eq.1) then
C------- Translate R- or SBEND to BEND
        if(ik.eq.2) then
C--------- For RBEN -> BEND, uncomment the following line "goto 21", otherwise will be RBEN -> MULTIPOL
cc         goto 21
        elseif(ik.eq.3) then
C--------- SBEN -> BEND
          goto 3
        endif
      endif
C----- R- or SBEN -> MULTIPOL
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif
      ss = ss+xl
      xlarc = xl
      if(ang.ne.0.d0) then
        ro = xlarc /ang
        xlmag = 2.d0 * ro * sin (ang/2.d0)
        b1 = bro / ro                              ! dipole field (T)
      else
        xlmag=xl
        b1=0.d0
      endif
      
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
      endif
! quadrupole field at .1 m (T)
      b2 = ak1*bro * x10/100.                     
! sextupole field at .1 m (T)
      b3 = ak2 * bro *(x10/100.)**2 / 2.d0
C         write(*,*) ak2, bro, b3, x10

      if(b1*b1+b2*b2+b3*b3.eq.0.d0) goto 1

      write(lw,fmt='(A,T12,A,T22,A)') '''MULTIPOL''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .Dip'
      write(lw,fmt='(F9.4,F7.2,3F13.8,7F3.0)')
     >xlmag*cm, x10, b1*10.d0, b2*10.d0, b3*10.d0, x0,x0,x0,x0,x0,x0,x0

      txt=txffd
      if(frf .eq. 0.) txt='0. 0. '//txt
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A,2X,A)') txtsD,name
      write(lw,*)  kpos,xce,yce,ale
C      if(kpos.eq.3) then
C        write(lw,fmt='(I,1P,3G18.10)') kpos,zero,zero,zero
C        write(lw,fmt='(I,A,1P,2G18.10)') kpos,' 0. ',-off*100.,zero
C        write(lw,fmt='(I,A,1P,2G18.10)') kpos,' 0. ',-off*100.,ang/2.
C      else
C        write(lw,fmt='(I,1P,3E16.8)') kpos,xce,yce,ale
C      endif
      prec = kley(1:4)
      goto 99    
 
 21   continue
C----- RBEN -> BEND
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif

      ss = ss+xl
      ro=xl/2./sin(ang/2.)
      b=ang * bro / xl *10.  !bro/ro*tkg
C      xxl=ro*ang*cm
      ro=ro*cm

      if(b.eq.0.d0) goto 1

      write(lw,fmt='(A,T12,A,T22,A)') '''BEND''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .Bend'
      write(lw,fmt='(2F14.7,F15.8)') xl*cm,tilt,b 

      te=ang/2.+e1
      if(frf .eq. 0.) then
        write(txt,*) x0,x0,te
      else
        write(txt,*) x20,x8,te
      endif
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      ts=ang/2.+e2
      if(frf .eq. 0.) then
        write(txt,*) x0,x0,ts
      else
        write(txt,*) x20,x8,ts
      endif
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
C      step = ro * sin(ang/2.) * 2.
C      write(lw,fmt='(F12.6,A4)') step / 10.,'  Bend'
      write(lw,fmt='(A,2X,A)') ' #20|60|20    Rbend',name
      if(kpos.eq.3) then
        write(lw,fmt='(A)') '3 0. 0. 0.'
      else
        write(lw,fmt='(I,1P,3G18.10)') kpos,xce,yce,ale
      endif
      prec = kley(1:4)
      goto 99

 3    continue
C----- SBEN -> BEND
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif

      ss = ss+xl
      ro=xl/ang
      xxl = 2.d0*ro*sin(ang/2.d0)
           xxl=xl
      b=bro/ro*tkg

      if(b.eq.0.d0) goto 1

      write(lw,fmt='(A,T12,A,T22,A)') '''BEND''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .Bend'
      write(lw,fmt='(2F14.7,F10.4)') xxl*cm,tilt,b 
      te=e1
      if(frf .eq. 0.) then
        write(txt,*) x0,x0,te
      else
        write(txt,*) x20,x8,te
      endif
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
      ts=e2
      if(frf .eq. 0.) then
        write(txt,*) x0,x0,ts
      else
        write(txt,*) x20,x8,ts
      endif
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
      write(lw,fmt='(A,2X,A)') ' #20|60|20    Sbend',name
      if(kpos.eq.3) then
        write(lw,fmt='(A)') '3 0. 0. 0.'
      else
        write(lw,fmt='(I,1P,3G18.10)') kpos,xce,yce,ale
      endif
      prec = kley(1:4)
      goto 99

 4    continue
C----- QUAD
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif
      ss = ss+xl

      if(ak1.eq.0.d0) goto 1

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
      write(lw,fmt='(A,2X,A)') txtsQ,name
      write(lw,fmt='(A)') '1 0. 0. 0.'
      prec = kley(1:4)
      goto 99
 5    continue
C----- SEXT ***
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif

      ss = ss+xl

      b3 = 10.d0 * ak2 * bro *(x10/100.)**2 / 2.d0
       b3=b3 * fb3

      if(b3.eq.0.d0) goto 1

      write(lw,fmt='(A,T12,A,T22,A)') '''MULTIPOL''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .Sext'
      write(lw,fmt='(F10.4,F6.2,2F10.4,F16.10,7F4.1)')
     >xl*cm,x10,x0,x0,b3,x0,x0,x0,x0,x0,x0,x0
      write(lw,fmt='(A)') '0.00 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') '0.00 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A,2X,A)') '#20|40|20  Sext',name
      write(lw,fmt='(A)') '1 0. 0. 0.'
      prec = kley(1:4)
      goto 99
 6    continue
C----- OCTU ***
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif
      ss = ss+xl

      if(b4.eq.0.d0) goto 1

      write(lw,fmt='(A,T12,A,T22,A)') '''MULTIPOL''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .Octu'
      write(lw,fmt='(F10.4,F6.2,3F10.4,F16.10,6F4.1)')
     >xl*cm,x10,x0,x0,x0,b4,x0,x0,x0,x0,x0,x0
      write(lw,fmt='(A)') '0.000 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') '0.000 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A,2X,A)') '#20|40|20  Octu',name
      write(lw,fmt='(A)') '1 0. 0. 0.'
      prec = kley(1:4)
      goto 99
 7    continue
C----- MULT
      if(xl .eq. 0) then
        it = it - 1
        goto 99
      endif
      ss = ss+xl

      if(ak1.eq.0.d0) goto 1

      write(lw,fmt='(A,T12,A,T22,A)') '''MULTIPOL''',kley,name
      write(lw,fmt='(I1,A)') i0,'  .Mult'
      write(lw,fmt='(F10.4,F6.2,2F16.10,8F4.1)')
     >xl*cm,x10,x0,ak1*bro,x0,x0,x0,x0
      write(lw,fmt='(A)') '00.00 00.0  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') '00.00 00.0  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A,2X,A)') '#20|60|20  Mult',name
      write(lw,fmt='(F8.2)') ten
      write(lw,fmt='(A)') '1 0. 0. 0.'
      prec = kley(1:4)
      goto 99
 8    continue
C----- SROT ***
      write(lw,fmt='(A,T12,A,T22,A)') '''TRAROT''',kley,name
      prec = kley(1:4)
      goto 99
 9    continue
C----- YROT ***
      write(lw,fmt='(A,T12,A,T22,A)') '''TRAROT''',kley,name
      prec = kley(1:4)
      goto 99
 10   continue
C----- MARK ***
c      it = it - 1
      write(lw,fmt='(A,T12,A,T22,A)') '''MARKER''',kley,name
      prec = kley(1:4)
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
      b = b*cik
      ss = ss+xl

      if(b.eq.0.d0) goto 1

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
      write(lw,fmt='(A)') '#20|20|20  Kick'
      write(lw,fmt='(A)') '1 0. 0. 0.'
      prec = kley(1:4)
      goto 99
 12   continue
C----- HMON, VMON
c      if(xl .eq. 0.d0) then
c        it = it - 1
c        goto 99
c      else
c        ss = ss+xl
c        write(lw,fmt='(A,T12,A,T22,A)') '''DRIFT''',kley,name
c        write(lw,fmt='(F12.4)') xl*cm
c        prec = kley(1:4)
c      endif
c      goto 99
      goto 1
 13   continue
C----- RCOL, ECOL
c      if(xl .eq. 0.d0) then
c        it = it - 1
c        goto 99
c      else
c        ss = ss+xl
c        write(lw,fmt='(A,T12,A,T22,A)') '''DRIFT''',kley,name
c        write(lw,fmt='(F12.4)') xl*cm
c        prec = kley(1:4)
c      endif
c      goto 99
      goto 1

 99   continue
       write(6,fmt='(A,T12,A,)') kley,name
       return

      entry lmntff(kff)

      if(kff .eq.1) then
!-------- LHC
        txffd = ffdlhc
        txfd = fdlhc
        txffq = ffqlhc
        txfq = fqlhc        
        txtsD = sDlhc
        txtsQ = sQlhc
      else if(kff .eq.2) then
!-------- Recycler      
        txffd = ffdrec
        txfd = fdrec
        txffq = ffqrec
        txfq = fqrec
        txtsD = stpD
        txtsQ = stpQ
      else if(kff .eq.3) then
!-------- muon collider 
        txffq = ffqmua
        txfq = fqmua
        txffd = ffdlhc
        txfd = fdlhc
        txtsD = stpD
        txtsQ = stpQ
      else
        txffd = ffd
        txfd = fd
        txffq = ffq
        txfq = fq
        txtsD = stpD
        txtsQ = stpQ
      endif

c      write(*,fmt='(A)') txffd ,txfd,txffq ,txfq

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

      entry lmnkbm(kbmi)
      kbm =kbmi 
      return

      entry lmnkic(ciki)
      cik=ciki
      return

      entry lmncon(ny)
      condri=ny.eq.'y'
      return

      entry lmnfb3(fb3i)
      fb3=fb3i
      return

      end
