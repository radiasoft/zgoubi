C23456789012345678901234567890123456789012345678901234567890123456789012
      implicit double precision (a-h,o-z)
      character fnr*80, fnw*80
      character txt500*500,keyword*20,txt20*20
      character txt132*132,txt14*14,txt40*40
      parameter(lr=9, lw=10)
      parameter (mxl=15000)
      dimension iq(mxl)
      character kley*4, name*20
      parameter(mxk=19)
      character kle(mxk)*4, ny*1
      character*35 warn
      logical TOMANY
      DIMENSION BK1(6), BK2(6)

      logical ok, strcon
      integer debstr, finstr

      data kle/ 'DRIF', 'RBEN', 'SBEN', 'QUAD', 'SEXT', 'OCTU', 'MULT',
     >'SROT', 'YROT', 'MARK', 'KICK', 'HMON', 'VMON','HKIC', 'VKIC', 
     >'SOLE', 'RCOL', 'MATR', 'MONI' /

C      write(*,*) ' name of the twiss file to be translated:'
C      read(*,fmt='(A)') fnr
C      write(*,*) ' name of the zgoubi data file:'
C      read(*,fmt='(A)') fnw

C fnr is a 'twiss' type of file from MAD. 
Case MAD8 : as resulting from a "twiss, save"
Case MADX : as resulting from a "select" w appropriate argument list
      fnr='madzg.in'
      if(fnr.ne.'madzg.in') then
        write(*,*)
        write(*,*) 'fnr MUST be  ''madzg.in'' '
        write(*,*)
      endif
      fnw='trad.out'

      OPEN(UNIT=lr,FILE=fnr,status='old',err=97)
 80   OPEN(UNIT=lw,FILE=fnw,STATUS='NEW',IOSTAT=IOS)
      IF(IOS .NE. 0) THEN
        CLOSE(lw,STATUS='DELETE')
        OPEN(UNIT=lw,FILE=fnw,STATUS='OLD')
        CLOSE(lw,STATUS='DELETE')
        GOTO 80
      ENDIF
     
 8    write(*,*) 
      write(*,*) 
     >'Centered or long-shifted or short-shftd multipole (1/2/3)?'
C      read(*,*,err=8) kcs
          kcs = 2
      write(*,*) ' kcs = ', kcs
      if( kcs .lt. 1 .or. kcs .gt. 3) goto 8
      call lmnkcs(kcs)

      write(*,*) 
      write(*,*) ' Fringe field coeff :',
     >             ' lhc/recycler/musr/default (1/2/3/99)'
C 2    read(*,*,err=2) kff
         kff = 99
C      if( kff .ne. 1 .and. kff .ne. 2) kff = 3
      write(*,*) ' Fringe field coeff : option set to ',kff
      call lmntff(kff)

      write(*,*) 
      write(*,*) ' Fringe fields on/off/default (1/0/99) :'
c 3    read(*,*,err=3) frf
          frf  = 0
      if( frf .ne. 1.d0) frf = 0.d0
      if(frf.eq.1.d0) write(*,*) ' Fringe fields on '
      if(frf.eq.0.d0) write(*,*) ' Fringe fields off '

      write(*,*) 
      write(*,*) ' Option KPOS in BEND and MULTIPOLE - ', 
     > 'Option 4 is to change KPOS=3 into the equivalent CHANGREF',
     >                      '1/2/3/4/default (1/2/3/4/99) :'
c 4    read(*,*,err=4) kpos
             kpos = 3
      if( kpos .lt. 1 .or. kpos .gt. 4) kpos=2
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
      write(*,*) ' Introduce CHANGREF at main bends (n/y) :'
      read(*,*,err=51) ny
 51   if(ny.ne.'y') ny='n' 
      write(*,*) '               option set to ',ny
      call lmndr8(ny)

      write(*,*)
      write(*,*) ' Convert DRIFT to zero-field MULTIPOLE (n/y) :'
c 6    read(*,*,err=6) ny
              ny = 'n'
      if(ny.ne.'y') ny='n' 
      write(*,*) '               option set to ',ny
      call lmndri(ny)

      write(*,*)
      write(*,*) ' Translate S/RBEND to MULTIPOL / AGSMM KPOS=3 / AGSMM'
     > ,' KPOS=4 (2 / 3 / 4) :'
      read(*,*,err=7) kbm
C            kbm = 3
 7    continue
C      if( kbm .ne. 2   ) kbm = 1
      if(kbm.eq.1   ) write(*,*) ' BEND will be used'
      if(kbm.eq.2) write(*,*) ' MULTIPOLE will be used for S/RBEND'
      if(kbm.eq.3) write(*,*) ' AGSMM w/ KPOS=3 will be used '
      if(kbm.eq.4) write(*,*) ' AGSMM w/ KPOS=4 will be used '
      call lmnkbm(kbm)

      write(*,*)
      write(*,*) ' Translate QUAD to AGSQUAD  / MULTIPOL (1/2)'
      write(*,*) ' WARNING : AGSQUAD will have currents set to zero.'
      read(*,*,err=9) kquad
c      kquad = 2
 9    continue
      if( kquad .ne. 2   ) kquad = 1
      if(kquad.eq.1) write(*,*) ' AGSQUAD will be used '
      if(kquad.eq.2) write(*,*) ' MULTIPOLE will be used '
      call lmnkq(kquad)


C Get momentum
 32   continue
        read(lr,fmt='(A)',end=85,err=99) txt500
        ok =       strcon(txt500,'@ MASS',
     >                                      IS) 
      if(.not. ok) goto 32        

      read(txt500(finstr(txt500)-17:finstr(txt500)),*) amass     ! GeV/c2
        read(lr,fmt='(A)',end=85,err=99) txt500
      read(txt500(finstr(txt500)-17:finstr(txt500)),*) nq
        read(lr,fmt='(A)',end=85,err=99) txt500
      read(txt500(finstr(txt500)-17:finstr(txt500)),*) energ   ! GeV
        read(lr,fmt='(A)',end=85,err=99) txt500
      read(txt500(finstr(txt500)-17:finstr(txt500)),*) pmom    ! GeV/c

c       write(*,*) pmom
       
      it = 0
      call objet(lw,it,pmom*1.d9,nq,amass*1d9,
     >                           bro)        
      bro = 1.d0 

      CALL strng(pmom/dble(nq),
     >                        bK1,bK2)

        write(*,*) ' bro, p : ',bro, pmom
        write(89,*) ' madzg K1 : ',bk1
        write(89,*) ' madzg K2 : ',bk2


C Swallow rest of top of madzg.in file
 33   continue
        read(lr,fmt='(A)',end=85,err=99) txt500
C        write(*,*) txt500
        ok =       strcon(txt500,'KEYWORD',
     >                                    IS) 
      if(.not. ok) goto 33
      read(lr,fmt='(A)',end=85,err=99) txt500
c      write(*,*) txt500

      write(*,*) 
      write(*,*) ' Now translating. Busy...'

      ir = 0
      noel = 0
      sxl = 0.d0
      s1 = 0.d0
      s = 0.d0
 86   continue
      
        read(lr,fmt='(a)',end=85,err=99) txt500
C        write(*,*) txt500
        read(txt500,fmt='(24e19.0,a)',end=85,err=99)
     >  xl, ang, ak1,ak2,ak3,ak4, tilt, e1,e2, h1,h2,sCntr,
     >  alfx,betx, alfy,bety, xco,yco, Dx,Dxp, Dy,Dyp, xmu,ymu, txt40

        read(txt40,*,end=85,err=99) txt20,keyword
        name = txt20
        kley = keyword
c        write(*,*) '-----------------------------------'
c        write(*,*) kley, '  ', name, xl, ang, ak1, ak2

        s = s+xl

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
     >    (lw,bro,frf,ik,noel,name,kley,xl,ang,ak1,ak2,ak3,ak4,
     >                             bk1,bk2,   tilt,e1,e2,h1,h2,
     >                                      it)
        ds = s - s1
        if(abs(ds-xl).gt.1d-3) then
          warn = 'WARNING : ds .ne. xl ; xl set to ds'
          xl = ds
          name = name(debstr(name):finstr(name))//' *** '
          kley = 'DRIF'
          it = it-1
          call lmnt
     >    (lw,bro,frf,ik,noel,name,kley,xl,ang,ak1,ak2,ak3,ak4,
     >                             bk1,bk2,     tilt,e1,e2,h1,h2,
     >                                      it)
        else
          warn = ' '
        endif
        s1 = s
        sxl = sxl + xl
        write(88,fmt='(1p,5(e16.8,1x),2(A16,1x),A35)') 
     >                         s,sxl,s-sxl,xl,ds,name,kley,warn

      goto 86

 85   continue
      write(lw,fmt='(A,T111,I6)') '''MATRIX''',it
      write(lw,fmt='(A)')   '2  11 '
      write(lw,fmt='(A,T111,I6)') '''TWISS''',it
      write(lw,fmt='(A)')   '1  1.  1. '
c      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''END''',it
C      it = it + 1
      write(lw,fmt='(///,A)') '''REBELOTE'''
      write(lw,fmt='(A)')   '600 0.1 99'
      write(*,*) ' Read ',ir,' elements from the madzg.in file ',fnr
      write(*,*) ' Translated to',it,
     >' elements into the zgoubi data file "trad.out"', fnw
      write(*,*) ' end of madzg.in file '
      goto 999

 97   write(*,*) ' error open madzg.in file'
      goto 999
 99   write(*,*) ' error during read of madzg.in file'
      goto 999

 999  continue
      close(lr)
      close(lw)
      stop
      end
        
      subroutine objet(lw,it,p,nq,am,
     >                            bro)
      implicit double precision (a-h,o-z)
     
 1    continue

        bro = p/dble(nq) /2.99792458d8    ! T.m 
C        am = 938.27203d6
        gamma = sqrt(p*p + am*am)/am
      write(*,*) ' Problem rigidity (T.m), mass (GeV), G.gamma : ',
     >   bro, am/1.d9, 1.79284735d0*gamma
          
      write(lw,*) 'Generated by MADX -> Zgoubi translator'
      write(lw,fmt='(A,T111,I6)') '''OBJET''',it
      write(lw,FMT='(F17.4,A,f17.4)') bro*1.d3,
     >'     reference rigidity (kG.cm),  G.gamma =  ',1.79284735d0*gamma
      write(lw,fmt='(A)') '5' 
C      write(lw,fmt='(A)') '.001 .001 .001 .001 0. .0001  '
      write(lw,fmt='(A)') '.01 .01 .01 .01 0. .001  '
      write(lw,fmt='(A)') '0. 0. 0. 0. 0. 1.'

      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''PARTICUL''',it
      write(lw,FMT='(1p,2(e15.7,1x),A)') am/1.d6, 
     >      dble(nq)*1.602176487E-19,' 1.7928474 0 0' 

      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''FAISCEAU''',it

c      it = it + 1
c      write(lw,fmt='(A,T111,I6)') '''OPTICS''',it
c      write(lw,FMT='(A)') 
c     > ' 2   Print out transport coeffs to zgoubi_MATRIX_out' 

      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''SCALING''' 
      write(lw,fmt='(A)') '1  1 '
c      write(lw,fmt='(A)') 'AGSMM' 
c      write(lw,fmt='(A)') '   -1   ' 
c      write(lw,fmt='(A)') '   1.   ' 
C      write(lw,fmt='(1P,E16.8)') BRO
c      write(lw,fmt='(A)') '1  '
C      write(lw,fmt='(A)') 'MULTIPOL   QH_*'
C      write(lw,fmt='(A)') '-1         '
C      write(lw,fmt='(1P,E16.8)') BRO
C      write(lw,fmt='(A)') '1      '
C      write(lw,fmt='(A)') 'MULTIPOL   QV_*'
C      write(lw,fmt='(A)') '-1         '
C      write(lw,fmt='(1P,E16.8)') BRO
C      write(lw,fmt='(A)') '1      '
C      write(lw,fmt='(A)') 'MULTIPOL   QJUMI QJUMJ'
C      write(lw,fmt='(A)') '0  -87  '
C      write(lw,fmt='(1P,E16.8)') BRO
C      write(lw,fmt='(2A)') '43 .237 40    Start Gg, dN, #Turns/ramp'
      write(lw,fmt='(A)') 'MULTIPOL '
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8)') BRO
      write(lw,fmt='(A)') '1      '

      return
      end 

      subroutine lmnt
     >(lw,bro,frf,ik,noel,name,kley,xl,ang,ak1,ak2,ak3,ak4,
     >                       bk1,bk2,   tilt,e1,e2,h1,h2,
     >                                it)
      implicit double precision (a-h,o-z)
      dimension bk1(6),bk2(6)
      character kley*(*)

      character(*) name
      character ny*1
      
      character*80 txfd, txfq,fqlhc , fdlhc ,fqrec , fdrec,fd,fq
      character*80 txffd,txffq,ffqlhc,ffdlhc,ffqrec,ffdrec,ffd,ffq
      character*80 ffdmu, fdmu, ffqmua,ffqmus,fqmua
      save txfd,txfq,txffd, txffq
      logical drimul, chgrf
      save drimul, chgrf

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

      parameter (cm=100.d0, t2kg=10.d0)
      parameter(i0=0, i1=1, i2=2, i4=4, i6=6)
      parameter(x0=0.d0,x1=1.d0,x2=2.d0,x3=3.d0,x4=4.d0,x5=5.d0,x6=6.d0)
      parameter(x7=7.d0,x8=8.d0,x10=10.d0,x20=20.d0,x999=999.d0)

      character txt*80, txtfrm*23
      save kpos, kbm, kcs, kquad

      integer finstr

      pi = 4.d0*atan(1.d0)
      zero=0.d0
      txtfrm = '(A,T12,A,T24,A,T111,i6)'

C----- 1    2    3    4    5    6    7    8    9    10   11   12   13
C----- DRIF RBEN SBEN QUAD SEXT OCTU MULT SROT YROT MARK KICK HMON VMON
C----- 14   15   16   17   18   19
C----- HKIC VKIC SOLE RCOL MATR MONI
       goto (
     > 1,   2 ,  2 ,  4,   5,   6,   7,   8,   9,   10,  11,  12,  12,
     > 11,  11,  13,  1,  18,   12) ik

 1    continue
C----- DRIF      
      if(xl .eq. 0.d0) then
        write(lw,fmt=txtfrm)'''MARKER''',name,kley,it
        goto 99
      else
        if(drimul) then
          dum=1.d-20
          write(lw,fmt=txtfrm) '''MULTIPOL''',name,kley,it
          write(lw,fmt='(I1,A)') i0,'  .Drift'
          write(lw,fmt='(1P,E13.6,E11.4,E12.4,9E8.1)')
     >    xl*cm,x10,x0,dum,x0,x0,x0,x0,x0,x0,x0,x0
          txt = '0. 0. '//txffq
          write(lw,fmt='(A)') txt
          write(lw,fmt='(A)') txfq
          write(lw,fmt='(A)') txt
          write(lw,fmt='(A)') txfq
          write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
          write(lw,fmt='(A)') nint(xl*cm), '  Drift'
          write(lw,fmt='(A)') '1 0. 0. 0.'
        else
          write(lw,fmt=txtfrm) '''DRIFT''',name,kley,it
          write(lw,fmt='(F12.4)') xl*cm
        endif
      endif
      goto 99

 2    continue
       if(kbm.eq.1) then
C-------- Translate R- or SBEND to BEND
         if(ik.eq.2) then
C---------- RBEN -> BEND
           goto 21
         elseif(ik.eq.3) then
C---------- SBEN -> BEND
           goto 3
         endif
       elseif(kbm.eq.2) then
C----- R- or SBEN -> MULTIPOL 
         if(xl .eq. 0.d0) then
           write(lw,fmt=txtfrm) '''MARKER''',name,kley,it
           goto 99
         endif
       elseif(kbm.ge.3) then
C----- R- or SBEN -> AGSMM
         goto 22
       endif

      if(ik.eq.2) then
C RBEN
C MADX
        xlarc = xl
c MAD8        xlarc = xl/(2.d0*sin(ang/2.d0)) * ang
      elseif(ik.eq.3) then
C SBEN
        xlarc = xl
      else
        stop ' sbr lmnt, no such option ik '
      endif
      if(ang.ne.0.d0) then
        ro = xlarc /ang
        xlmag = 2.d0 * ro * sin (ang/2.d0)
        b1 = bro / ro                              ! dipole field (T)
        if(kpos.eq.3) then
C          ale = -0.5d0 * ang
          ale = 0.5d0 * ang
          xce=0.d0
          yce=0.d0
        elseif(kpos.eq.4) then
         ale = -0.5d0 * ang * 180.d0/pi
         yce = cm*ro*(1-cos(ale*pi/180.d0))
        endif
      else
        xlmag=xl
        b1=0.d0
          ale=0.d0
          xce=0.d0
          yce=0.d0
      endif

      if(kpos.eq.4) then
       kley = ' '
       name = ' '
       it = it+1
       write(lw,fmt=txtfrm) '''CHANGREF''',name,kley,it
       write(lw,fmt='(F13.9,3X,F13.9,3X,F13.9)') 0.d0,0.d0,ale 
       it = it+1
       write(lw,fmt=txtfrm) '''CHANGREF''',name,kley,it
       write(lw,fmt='(F13.9,3X,F13.9,3X,F13.9)') 0.d0,yce,0.d0
      endif
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .Dip'

! quadrupole field (T) at .1 m , and shift
        last = finstr(name)
        if(kcs .eq. 1) then
C centered dipole model
          sag = -0.  !21906d0 * 2.54d0   /16.d0*18.d0
          sagS = -0.  !15472d0 * 2.54d0  /16.d0*18.d0
          sgn = -1.d0
          fb1 = 1.d0   !/1.0028d0
          fb2 = 1.d0/fb1
          if    (name(last-1:last) .eq. 'BF') then
            sag = sagS
            xlmag = 79.d0 * 0.0254d0
            b3 = bk2(1) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(1) 
            b2 = bk1(1) - bk2(1)*(sag/2.d0)/100.d0
          elseif(name(last-1:last) .eq. 'CD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(2) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(2) 
            b2 = bk1(2) - bk2(2)*(sag/2.d0)/100.d0 
          elseif(name(last-1:last) .eq. 'AF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(3) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(3) 
            b2 = bk1(3) - bk2(3)*(sag/2.d0)/100.d0 
          elseif(name(last-1:last) .eq. 'BD') then  
            sag = sagS
            xlmag = 79.d0 * 0.0254d0
            b3 = bk2(4) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(4) 
            b2 = bk1(4) - bk2(4)*(sag/2.d0)/100.d0 
          elseif(name(last-1:last) .eq. 'CF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(5) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(5) 
            b2 = bk1(5) - bk2(5)*(sag/2.d0)/100.d0 
          elseif(name(last-1:last) .eq. 'AD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(6) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(6) 
            b2 = bk1(6) - bk2(6)*(sag/2.d0)/100.d0 
          endif
        elseif(kcs .eq. 2) then
C long-shifted dipole model
          fyce = 1.d0
          sgn = 1.d0
          b1 = 0.d0
          fb1 = 1.d0
          fb2 = 1.d0
          if    (name(last-1:last) .eq. 'BF') then
            xlmag = 79.d0 * 0.0254d0
            b3 = -bk2(1) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (23.17d0 )
            b2 = bk1(1) - bk2(1)*yce/100.d0
          elseif(name(last-1:last) .eq. 'CD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = -bk2(2) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-24.07d0 )
            b2 = bk1(2) - bk2(2)*yce/100.d0 
          elseif(name(last-1:last) .eq. 'AF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = -bk2(3) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (23.d0 )
            b2 = bk1(3) - bk2(3)*yce/100.d0 
          elseif(name(last-1:last) .eq. 'BD') then  
            xlmag = 79.d0 * 0.0254d0
            b3 = -bk2(4) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-23.93d0 )
            b2 = bk1(4) - bk2(4)*yce/100.d0 
          elseif(name(last-1:last) .eq. 'CF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = -bk2(5) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (23.d0 )
            b2 = bk1(5) - bk2(5)*yce/100.d0 
          elseif(name(last-1:last) .eq. 'AD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = -bk2(6) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-24.07d0 )
            b2 = bk1(6) - bk2(6)*yce/100.d0 
          endif
c          fyce = 1.d0
c          sgn = 1.d0
c          fb1 = 1.d0
c          fb2 = 1.d0
c          if    (name(last-1:last) .eq. 'BF') then
c            xlmag = 79.d0 * 0.0254d0
c            b3 = -bk2(1) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (23.17d0 )
c            b2 = bk1(1) 
c          elseif(name(last-1:last) .eq. 'CD') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(2) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (-24.07d0 )
c            b2 = bk1(2) 
c          elseif(name(last-1:last) .eq. 'AF') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(3) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (23.d0 )
c            b2 = bk1(3) 
c          elseif(name(last-1:last) .eq. 'BD') then  
c            xlmag = 79.d0 * 0.0254d0
c            b3 = -bk2(4) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (-23.93d0 )
c            b2 = bk1(4) 
c          elseif(name(last-1:last) .eq. 'CF') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(5) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (23.d0 )
c            b2 = bk1(5) 
c          elseif(name(last-1:last) .eq. 'AD') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(6) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (-24.07d0 )
c            b2 = bk1(6) 
c          endif
c          b1 = b1 + b2*yce/100.d0 - b3 * yce*yce/1.d4
c          yce = -yce
        elseif(kcs .eq. 3) then
C short-shifted dipole model
          fyce = 1.d0
          sgn = -1.d0
          fb1 = 1.d0
          fb2 = 1.d0
          if    (name(last-1:last) .eq. 'BF') then
c               write(*,*) xl,xlmag,xlarc
            xlmag = 79.d0 * 0.0254d0
            b3 = bk2(1) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.3942767240d0)
            b2 = bk1(1) 
          elseif(name(last-1:last) .eq. 'CD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(2) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.5539466648d0)
            b2 = bk1(2) 
          elseif(name(last-1:last) .eq. 'AF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(3) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.5589609130d0)
            b2 = bk1(3) 
          elseif(name(last-1:last) .eq. 'BD') then  
            xlmag = 79.d0 * 0.0254d0
            b3 = bk2(4) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.3917462578d0)
            b2 = bk1(4) 
          elseif(name(last-1:last) .eq. 'CF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(5) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.5589406517d0)
            b2 = bk1(5) 
          elseif(name(last-1:last) .eq. 'AD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(6) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.5539271156d0)
            b2 = bk1(6) 
          endif
        else
          write(*,*) ' No such possibility kcs = ',kcs
          stop
        endif

! sextupole field at .1 m (T)
C      b3 = - ak2 * bro *(x10/100.d0)**2 / 2.d0

c      write(*,*) xlmag*cm, x10, x0, b2, 2.*b3*10.d0

      write(lw,fmt='(E16.8,1X,F8.4,3(1X,F14.10),1X,7F4.1)')
     >  xlmag*cm,x10,
C     >  sgn*fb1*b1*10.d0, fb2*b2, b3*10.d0*sgn, x0,x0,x0,x0,x0,x0,x0
     >  -sgn*fb1*b1*10.d0,fb2*b2,b3*10.d0*sgn,x0,x0,x0,x0,x0,x0,x0

      txt=txffd
      if(frf .eq. 0.) txt='0. 0. '//txt
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
        write(lw,fmt='(2A)') ' 1.  Dip ',name
      if(kpos.ne.4) then
        ale =  -ale
C        write(*,*)  kpos,xce,yce,ale
        write(lw,fmt='(i1,2(2x,F11.7),2X,F13.9,a,f12.8,a)') 
C     >  kpos,xce,yce,sgn*ale,'  angle = ',2.*sgn*ale,' rad'
     >  kpos,xce,yce,-sgn*ale,'  angle = ',2.*sgn*ale,' rad'
      else
        write(lw,*)  1,0.d0,0.d0,0.d0
      endif

      if(kpos.eq.4) then
       name = ' '
       kley = ' '
       it = it+1
       write(lw,fmt=txtfrm) '''CHANGREF''',name,kley,it
       write(lw,fmt='(F13.9,3X,F13.9,3X,F13.9)') 0.d0,-yce,0.d0
       it = it+1
       write(lw,fmt=txtfrm) '''CHANGREF''',name,kley,it
       write(lw,fmt='(F13.9,3X,F13.9,3X,F13.9)') 0.d0,0.d0,ale 
      endif

      goto 99
 
 22   continue
C----- R- or SBEN -> AGSMM
      if(ik.eq.2) then
C RBEN
C MADX
        xlarc = xl
c MAD8        xlarc = xl/(2.d0*sin(ang/2.d0)) * ang
      elseif(ik.eq.3) then
C SBEN
        xlarc = xl
      else
        stop ' sbr lmnt, no such option ik '
      endif
      if(ang.ne.0.d0) then
        ro = xlarc /ang
        xlmag = 2.d0 * ro * sin (ang/2.d0)
        b1 = bro / ro                              ! dipole field (T)
        if(kpos.eq.3) then
C          ale = -0.5d0 * ang
          ale = 0.5d0 * ang
          xce=0.d0
          yce=0.d0
        elseif(kpos.eq.4) then
         ale = -0.5d0 * ang * 180.d0/pi
         yce = cm*ro*(1-cos(ale*pi/180.d0))
        endif
      else
        xlmag=xl
        b1=0.d0
          ale=0.d0
          xce=0.d0
          yce=0.d0
      endif

c      if(kpos.eq.4) then
c       name = ' '
c       kley = ' '
c       it = it+1
cC       write(lw,fmt=txtfrm) '''CHANGREF''',kley,name,it
c       write(lw,fmt=txtfrm) '''CHANGREF''',name,kley,it
c       write(lw,fmt='(F13.9,3X,F13.9,3X,F13.9)') 0.d0,0.d0,ale 
c       it = it+1
cC       write(lw,fmt=txtfrm) '''CHANGREF''',kley,name,it
c       write(lw,fmt=txtfrm) '''CHANGREF''',name,kley,it
c       write(lw,fmt='(F13.9,3X,F13.9,3X,F13.9)') 0.d0,yce,0.d0
c      endif
      if(chgrf) then
       name = ' AGSMM_UP'
       kley = ' '
       it = it+1
       write(lw,fmt=txtfrm) '''CHANGREF''',name,it
       write(lw,fmt='(F13.9,3X,F13.9,3X,F13.9)') x0, x0, x0
      endif
C      write(lw,fmt=txtfrm) '''AGSMM''',kley,name,it
      write(lw,fmt=txtfrm) '''AGSMM''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'      .AGSDip'

! quadrupole field (T) at .1 m , and shift
        last = finstr(name)
        if(kcs .eq. 1) then
C centered dipole model
          sag = -0.  !21906d0 * 2.54d0   /16.d0*18.d0
          sagS = -0.  !15472d0 * 2.54d0  /16.d0*18.d0
          sgn = -1.d0
          fb1 = 1.d0   !/1.0028d0
          fb2 = 1.d0/fb1
          if    (name(last-1:last) .eq. 'BF') then
            sag = sagS
            xlmag = 79.d0 * 0.0254d0
            b3 = bk2(1) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(1) 
            b2 = bk1(1) - bk2(1)*(sag/2.d0)/100.d0
          elseif(name(last-1:last) .eq. 'CD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(2) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(2) 
            b2 = bk1(2) - bk2(2)*(sag/2.d0)/100.d0 
          elseif(name(last-1:last) .eq. 'AF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(3) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(3) 
            b2 = bk1(3) - bk2(3)*(sag/2.d0)/100.d0 
          elseif(name(last-1:last) .eq. 'BD') then  
            sag = sagS
            xlmag = 79.d0 * 0.0254d0
            b3 = bk2(4) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(4) 
            b2 = bk1(4) - bk2(4)*(sag/2.d0)/100.d0 
          elseif(name(last-1:last) .eq. 'CF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(5) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(5) 
            b2 = bk1(5) - bk2(5)*(sag/2.d0)/100.d0 
          elseif(name(last-1:last) .eq. 'AD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(6) *(x10/100.d0)**2 / 2.d0
            b1 = b1 - (sag/2.d0)/100.d0 * bk1(6) 
            b2 = bk1(6) - bk2(6)*(sag/2.d0)/100.d0 
          endif
        elseif(kcs .eq. 2) then
C long-shifted dipole model
          fyce = 1.d0
          sgn = 1.d0
          b1 = 0.d0
          fb1 = 1.d0
          fb2 = 1.d0
          if    (name(last-1:last) .eq. 'BF') then
            xlmag = 79.d0 * 0.0254d0
            b3 = -bk2(1) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (23.17d0 )
            b2 = bk1(1) - bk2(1)*yce/100.d0
          elseif(name(last-1:last) .eq. 'CD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = -bk2(2) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-24.07d0 )
            b2 = bk1(2) - bk2(2)*yce/100.d0 
          elseif(name(last-1:last) .eq. 'AF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = -bk2(3) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (23.d0 )
            b2 = bk1(3) - bk2(3)*yce/100.d0 
          elseif(name(last-1:last) .eq. 'BD') then  
            xlmag = 79.d0 * 0.0254d0
            b3 = -bk2(4) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-23.93d0 )
            b2 = bk1(4) - bk2(4)*yce/100.d0 
          elseif(name(last-1:last) .eq. 'CF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = -bk2(5) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (23.d0 )
            b2 = bk1(5) - bk2(5)*yce/100.d0 
          elseif(name(last-1:last) .eq. 'AD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = -bk2(6) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-24.07d0 )
            b2 = bk1(6) - bk2(6)*yce/100.d0 
          endif
c          fyce = 1.d0
c          sgn = 1.d0
c          fb1 = 1.d0
c          fb2 = 1.d0
c          if    (name(last-1:last) .eq. 'BF') then
c            xlmag = 79.d0 * 0.0254d0
c            b3 = -bk2(1) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (23.17d0 )
c            b2 = bk1(1) 
c          elseif(name(last-1:last) .eq. 'CD') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(2) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (-24.07d0 )
c            b2 = bk1(2) 
c          elseif(name(last-1:last) .eq. 'AF') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(3) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (23.d0 )
c            b2 = bk1(3) 
c          elseif(name(last-1:last) .eq. 'BD') then  
c            xlmag = 79.d0 * 0.0254d0
c            b3 = -bk2(4) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (-23.93d0 )
c            b2 = bk1(4) 
c          elseif(name(last-1:last) .eq. 'CF') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(5) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (23.d0 )
c            b2 = bk1(5) 
c          elseif(name(last-1:last) .eq. 'AD') then  
c            xlmag = 94.d0 * 0.0254d0
c            b3 = -bk2(6) *(x10/100.d0)**2 / 2.d0
c            yce = fyce * (-24.07d0 )
c            b2 = bk1(6) 
c          endif
c          b1 = b1 + b2*yce/100.d0 - b3 * yce*yce/1.d4
c          yce = -yce
        elseif(kcs .eq. 3) then
C short-shifted dipole model
          fyce = 1.d0
          sgn = -1.d0
          fb1 = 1.d0
          fb2 = 1.d0
          if    (name(last-1:last) .eq. 'BF') then
c               write(*,*) xl,xlmag,xlarc
            xlmag = 79.d0 * 0.0254d0
            b3 = bk2(1) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.3942767240d0)
            b2 = bk1(1) 
          elseif(name(last-1:last) .eq. 'CD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(2) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.5539466648d0)
            b2 = bk1(2) 
          elseif(name(last-1:last) .eq. 'AF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(3) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.5589609130d0)
            b2 = bk1(3) 
          elseif(name(last-1:last) .eq. 'BD') then  
            xlmag = 79.d0 * 0.0254d0
            b3 = bk2(4) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.3917462578d0)
            b2 = bk1(4) 
          elseif(name(last-1:last) .eq. 'CF') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(5) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.5589406517d0)
            b2 = bk1(5) 
          elseif(name(last-1:last) .eq. 'AD') then  
            xlmag = 94.d0 * 0.0254d0
            b3 = bk2(6) *(x10/100.d0)**2 / 2.d0
            yce = fyce * (-0.5539271156d0)
            b2 = bk1(6) 
          endif
        else
          write(*,*) ' No such possibility kcs = ',kcs
          stop
        endif

! sextupole field at .1 m (T)
C      b3 = - ak2 * bro *(x10/100.d0)**2 / 2.d0

c      write(*,*) xlmag*cm, x10, x0, b2, 2.*b3*10.d0

      write(lw,fmt='(A)') ' 3 0. 0. 0. 0. 0. '
      write(lw,fmt='(A)') ' 2 1 0. 1 0. '
C      write(lw,fmt='(E16.8,1X,F8.4,3(1X,F14.10),1X,7F4.1)')
C     >  xlmag*cm,x10,
C     >  sgn*fb1*b1*10.d0,fb2*b2,b3*10.d0*sgn,x0,x0,x0,x0,x0,x0,x0
C     >  -sgn*fb1*b1*10.d0,fb2*b2,b3*10.d0*sgn,x0,x0,x0,x0,x0,x0,x0

      txt=txffd
      if(frf .eq. 0.) txt='0. 0. '//txt
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
        write(lw,fmt='(2A)') ' 1.  AGSDip ',name
      if    (kbm.eq.3) then
        ale =  -ale
C        write(*,*)  kpos,xce,yce,ale
        write(lw,fmt='(i1,2(2x,F11.7),2X,F13.9,a,f12.8,a)') 
C     >  kpos,xce,yce,sgn*ale,'  angle = ',2.*sgn*ale,' rad'
     >  kpos,xce,yce,-sgn*ale,'  angle = ',2.*sgn*ale,' rad'
      elseif(kbm.eq.4) then
        write(lw,*) ' 4 0.  0.  0. 0. 0.'
      endif
      if(chgrf) then
       kley = ' '
       name = ' AGSMM_DW'
       it = it+1
       write(lw,fmt=txtfrm) '''CHANGREF''',name,kley,it
       write(lw,fmt='(F13.9,3X,F13.9,3X,F13.9)') x0, x0, x0
      endif

      if(kpos.eq.4) then
       kley = ' '
       name = ' '
       it = it+1
       write(lw,fmt=txtfrm) '''CHANGREF''',name,kley,it
       write(lw,fmt='(F13.9,3X,F13.9,3X,F13.9)') 0.d0,-yce,0.d0
       it = it+1
       write(lw,fmt=txtfrm) '''CHANGREF''',name,kley,it
       write(lw,fmt='(F13.9,3X,F13.9,3X,F13.9)') 0.d0,0.d0,ale 
      endif

      goto 99
 

 21   continue
C----- RBEN -> BEND
      if(xl * ang .eq. 0.d0) then
        if(xl .eq. 0.d0) then
          write(lw,fmt=txtfrm) '''MARKER''',name,kley,it
        else
          write(lw,fmt=txtfrm) '''DRIFT''',name,kley,it
          write(lw,fmt='(1P,E14.6)') xl*cm
        endif
        goto 99
      endif
c      if(xl .eq. 0.d0) then
c        write(lw,fmt=txtfrm) '''MARKER''',name,kley,it
c        goto 99
c      endif

      write(lw,fmt=txtfrm) '''BEND''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .Bend'
      ro=xl/2./sin(ang/2.)
      b=ang * bro / xl *10.  !bro/ro*t2kg
C      xxl=ro*ang*cm
      ro=ro*cm

C err. FM June 2009      write(lw,fmt='(2F14.7,F15.8)') xl,tilt,b 
      write(lw,fmt='(2F14.7,F15.8)') xl*cm,tilt,b 
      te=ang/2.+e1
      if(frf .eq. 0.) then
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,te
      else
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,te
      endif
      write(lw,fmt='(A)') txfd
      ts=ang/2.+e2
      if(frf .eq. 0.) then
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,ts
      else
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,ts
      endif
      write(lw,fmt='(A)') txfd
      step = ro * sin(ang/2.) * 2.d0 *cm
C      write(lw,fmt='(1P,E12.4,A)') dble(nint(step)),'  Bend'
      write(lw,fmt='(A,2X,A)') ' #10|10|10   Bend',name
      if(kpos.eq.3) then
        write(lw,fmt='(A)') '3 0. 0. 0.'
      else
        write(lw,fmt='(1P,I1,2X,3G18.10)') kpos,xce,yce,ale
      endif
      goto 99

 3    continue
C----- SBEN -> BEND
      if(xl * ang .eq. 0.d0) then
        if(xl .eq. 0.d0) then
          write(lw,fmt=txtfrm) '''MARKER''',name,kley,it
        else
          write(lw,fmt=txtfrm) '''DRIFT''',name,kley,it
          write(lw,fmt='(1P,E14.6)') xl*cm
        endif
        goto 99
      endif
c      if(xl .eq. 0.d0) then
c        write(lw,fmt=txtfrm) '''MARKER''',name,kley,it
c        goto 99
c      endif

      write(lw,fmt=txtfrm) '''BEND''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .Bend'
      ro=xl/ang
      xxl = 2.d0*ro*sin(ang/2.d0)
      b=bro/ro*t2kg
C              write(*,*) ' ********* tilt ',tilt
      write(lw,fmt='(1P,3(1X,E15.7))') xxl*cm,tilt,b 
      te=e1
      if(frf .eq. 0.) then
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,te
      else
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,te
      endif
      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
      ts=e2
      if(frf .eq. 0.) then
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,ts
      else
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,ts
      endif
      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
      step = ro * ang *cm
C      write(lw,fmt='(1P,E12.4,A)') dble(nint(step)),'  Bend'
      write(lw,fmt='(A,2X,A)') ' #10|10|10  Bend',name
      if(kpos.eq.3) then
        write(lw,fmt='(A)') '3 0. 0. 0.'
      else
        write(lw,fmt='(1P,I1,2X,3G18.10)') kpos,xce,yce,ale
      endif
      goto 99

 4    continue
C----- QUAD
      if(xl .eq. 0.d0) then
        write(lw,fmt=txtfrm) '''MARKER''',name,kley,it
        goto 99
      endif
      if   (kquad.eq.1) then
C-------- Translate to AGSQUAD
        write(lw,fmt=txtfrm) '''AGSQUAD''',name,kley,it
        write(lw,fmt='(I1,A)') i0,'  .Quad'  
        write(lw,fmt='(E14.6,1x,F8.4,1x,2(F16.10,1x),8(F4.1,1x))')
     >  xl*cm,x10,x0,x0,x0,x0,x0,x0
      elseif(kquad.eq.2) then
C-------- Translate to MULTIPOL
        write(lw,fmt=txtfrm) '''MULTIPOL''',name,kley,it
        write(lw,fmt='(I1,A)') i0,'  .Quad'  
        write(lw,fmt='(E14.6,F8.4,2F16.10,8F4.1)')
     >  xl*cm,x10,x0,ak1/xl*bro,x0,x0,x0,x0,x0,x0,x0,x0
      endif
      txt = txffq
      txt=txt//name
      if(frf .eq. 0.) txt = '0. 0. '//txt
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfq
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfq
      write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A,2X,A)') ' 1.0  Quad',name
      write(lw,fmt='(A)') '1 0. 0. 0.'
      goto 99

 5    continue
C----- SEXT ***
      if(xl .eq. 0.d0) then
        write(lw,fmt=txtfrm) '''MARKER''',name,kley,it
        goto 99
      endif
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .Sext'
      b3 = - 10.d0 * ak2/xl * bro *(x10/100.d0)**2 / 2.d0
C      write(*,*) ' SEXT         b3 = ',b3,' kG'
      write(lw,fmt='(E14.6,F8.4,2F10.4,F16.10,7F4.1)')
     >xl*cm,x10,x0,x0,b3,x0,x0,x0,x0,x0,x0,x0
      write(lw,fmt='(A)') '0.00 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') '0.00 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A,2X,A)') ' 1.00  Sext',name
      write(lw,fmt='(A)') '1 0. 0. 0.'
      goto 99

 6    continue
C----- OCTU ***
      if(xl .eq. 0.d0) then
        write(lw,fmt=txtfrm) '''MARKER''',name,kley,it
C        it = it - 1
        goto 99
      endif
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .Octu'
      write(lw,fmt='(E14.6,F8.4,3F10.4,F16.10,6F4.1)')
     >xl*cm,x10,x0,x0,x0,x0,x0,x0
      write(lw,fmt='(A)') '0.000 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') '0.000 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A,2X,A)') '1.000  Octu',name
      write(lw,fmt='(A)') '1 0. 0. 0.'
       goto 99

 7    continue
C----- MULT
      if(abs(ak1)+abs(ak2)+abs(ak3)+abs(ak4) .gt. 1.d-30) then
        xl = 1.d-5
      else
        write(lw,fmt=txtfrm)'''MARKER''',name,kley,it
C        it = it - 1
        goto 99
      endif
! quadrupole field at .1 m (kG)
      b2 = ak1 * bro / xl 
! sextupole field at .1 m (kG)
      b3 = ak2 * bro *(x10/100.d0)**2 / 2.d0  / xl * 10.d0
! octupole field at .1 m (kG)
      b4 = 0.d0
      b5 = 0.d0
      b6 = 0.d0
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .Mult'
      write(lw,fmt='(F12.6,F8.4,3F16.10,7F4.1)')
     >xl*cm,x10,x0,b2,b3,b4,b5,b6,x0,x0,x0,x0
c      write(*,*)
c     >b1,b2,b3,b4,b5,b6,bro,' b1,b2,b3,b4,b5,b6, bro'
      write(lw,fmt='(A)') '00.00 00.0  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') '00.00 00.0  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(I6,2(2X,A))') 2+nint(xl/10.d0),'  Mult',name
      write(lw,fmt='(A)') '1 0. 0. 0.'
       goto 99

 8    continue
C----- SROT ***
      write(lw,fmt=txtfrm) '''TRAROT''',name,kley,it
      goto 99

 9    continue
C----- YROT ***
      write(lw,fmt=txtfrm) '''TRAROT''',name,kley,it
      goto 99

 10   continue
C----- MARK ***
      write(lw,fmt=txtfrm) '''MARKER''',name,kley,it
C      it = it - 1
      goto 99

 11   continue
C----- KICK ***
      xxl = xl
      if(xl .eq. 0.d0) xxl = 1.e-6
      if(ik.eq.14) then
c--------- hkicker
        diptlt = 0.
        b = -bro*tilt/xxl*t2kg
      elseif(ik.eq.15) then
c--------- vkicker
        diptlt = pi/2.d0
        b = -bro*e1/xxl*t2kg
      endif
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .kicker'
      write(lw,fmt='(E14.6,F8.4,E14.6,9F4.1)')
     >xxl*cm,x10,b,x0,x0,x0,x0,x0,x0,x0,x0,x0
      txt = '.0 .0  1.00 0.00 0.00 0.00 0.00 0. 0. 0. 0.'
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(f12.9,A)') diptlt,' 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A)') '#10|10|10  Kick'
      write(lw,fmt='(A)') '1 0. 0. 0.'
      goto 99

 12   continue
C----- HMON, VMON
      write(lw,fmt=txtfrm) '''DRIFT''',name,kley,it
      write(lw,fmt='(F12.4)') xl*cm
      goto 99

 13   continue
C 'SOLENOID'                                                                            2
C   0
C   3000.  26.1   15.         ***************************         0.29866369046pourEps=.0223
C   50.  100.            100. 100. **************************
C #400|1200|400
C   1  0. 0. 0.
C----- SOLE -> SOLENOID
      if(xl .eq. 0.d0) then
        write(lw,fmt=txtfrm) '''MARKER''',name,kley,it
C        it = it - 1
        goto 99
      endif
      r0 = 111.11
      write(lw,fmt=txtfrm) '''SOLENOID''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .soleno'
      write(lw,fmt='(1P,G14.7,G11.4,G14.7)') xl*cm,r0,e1*bro*t2kg
      write(lw,fmt='(2F10.2)') r0, r0
      write(lw,fmt='(f10.4,2x,A,2X,A)') xl,'  Soleno',name
      write(lw,fmt='(A)') '1 0. 0. 0. '
      goto 99

 18   continue
C----- MATR
      if(xl .eq. 0.d0) then
        write(lw,fmt=txtfrm)'''MARKER''',name,kley,it
        goto 99
      else
        if    (name.eq.'CSNK2') then
          write(lw,fmt='(A)')   '''DRIFT''    DRIF      MATR CSNK'
          write(lw,fmt='(A)')   '10.0000' 
          write(lw,fmt='(A)')   '''MARKER''   MATR      CSNK '
          write(lw,fmt='(A)')   '''MAP2D''' 
          write(lw,fmt='(A)')   '0  0    '
          write(lw,fmt='(A)')   '.7e-3   100. 100. 100. '
          write(lw,fmt='(A)')   'HEADER_4 csnake   '
          write(lw,fmt='(A)')   '281 29 '
          write(lw,fmt='(A)')   'ref+sole.map2d' 
          write(lw,fmt='(A)')   '0 0 0 0   '
          write(lw,fmt='(A)')   '2     '
          write(lw,fmt='(A)')   '.1     '
C          write(lw,fmt='(A)')   '2  0.  .15  0.  0.'  
          write(lw,fmt='(A)')   '2  0.  2.5  0.  0.'  
          write(lw,fmt='(A)')   '''DRIFT''' 
          write(lw,fmt='(A)')   '10.0000   '  
        elseif(name.eq.'WSNK') then
          write(lw,fmt='(A)')  '''DRIFT''    DRIF      MATR WSNK  '
          write(lw,fmt='(A)')  '-50.                     '
          write(lw,fmt='(A)')  '''MARKER''   MATR      WSNK'
          write(lw,fmt='(A)')  '''MAP2D''                  '
          write(lw,fmt='(A)')  '0  0                     '
          write(lw,fmt='(A)')  '1.e1   100. 100. 100.    '
          write(lw,fmt='(A)')  'HEADER_4 wsnake          '
          write(lw,fmt='(A)')  '801 29                   '
          write(lw,fmt='(A)')  'table55.map2d            '
          write(lw,fmt='(A)')  '0 0 0 0                  '
          write(lw,fmt='(A)')  '2                        '
          write(lw,fmt='(A)')  '.1                       '
C          write(lw,fmt='(A)')  '2  0.  .15  0.  0.       '
          write(lw,fmt='(A)')  '2  0.  2.  0.  0.       '
          write(lw,fmt='(A)')  '''DRIFT''    DRIF      MATR WSNK'
          write(lw,fmt='(A)')  '-50.                     '
        endif

      endif
      goto 99

 99    continue
c       write(*,fmt='(A,T12,A)') name,kley
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

      entry lmndr8(ny)
      chgrf=ny.eq.'y'
      return

      entry lmnkbm(kbmi)
      kbm =kbmi 
      return

      entry lmnkq(kquadi)
      kquad =kquadi 
      return

      entry lmnkcs(kcsi)
      kcs =kcsi 
      return

      end
      FUNCTION DEBSTR(STRING)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DEBSTR
      CHARACTER * (*) STRING

C     --------------------------------------
C     RENVOIE DANS DEBSTR LE RANG DU
C     1-ER CHARACTER NON BLANC DE STRING,
C     OU BIEN 0 SI STRING EST VIDE ou BLANC.
C     --------------------------------------

      DEBSTR=0
      LENGTH=LEN(STRING)
1     CONTINUE
        DEBSTR=DEBSTR+1
        IF (STRING(DEBSTR:DEBSTR) .EQ. ' ') THEN
          IF(DEBSTR .EQ. LENGTH) THEN
            DEBSTR = 0
            RETURN
          ELSE
            GOTO 1
          ENDIF
        ENDIF

      RETURN
      END
      FUNCTION FINSTR(STRING)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER FINSTR
      CHARACTER * (*) STRING
C     --------------------------------------
C     RENVOIE DANS FINSTR LE RANG DU
C     DERNIER CHARACTER NON BLANC DE STRING,
C     OU BIEN 0 SI STRING EST VIDE ou BLANC.
C     --------------------------------------

      FINSTR=LEN(STRING)+1
1     CONTINUE
        FINSTR=FINSTR-1
        IF(FINSTR .EQ. 0) RETURN
        IF (STRING(FINSTR:FINSTR) .EQ. ' ') GOTO 1

      RETURN
      END
      subroutine strng(p,
     >                   aK1,aK2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension ak1(6), ak2(6)
      parameter (bdot=0.d0)

      conv = 0.0254d0
      ALA = 94d0*CONV        !length in meters
      ALB = 79d0*CONV
      ALC = 94d0*CONV

      UKAFM3 =  0.0D0
      UKAFM2 = -9.9628200D-6
      UKAFM1 = -2.7124400D-4
      UKAF0 =   0.0487709D0
      UKAF1 =   4.2519700D-5
      UKAF2 =  -1.5027310D-5 
      UKAF3 =   1.9900440D-6
      UKAF4 =  -1.2474890D-7
      UKAF5 =   3.7239700D-9
      UKAF6 =  -4.3303320D-11

      UKADM3 =  0.0D0
      UKADM2 =  1.0494200D-5
      UKADM1 =  2.7764700D-4
      UKAD0 =  -0.0487187D0
      UKAD1 =  -4.4223160D-5
      UKAD2 =   1.5432060D-5
      UKAD3 =  -2.0235660D-6
      UKAD4 =   1.2581210D-7
      UKAD5 =  -3.7265440D-9
      UKAD6 =   4.3053950D-11

      UKBFM3 =  0.0D0
      UKBFM2 = -1.0370300D-5
      UKBFM1 = -2.7908700D-4
      UKBF0 =   0.0485774D0
      UKBF1 =   4.1807340D-5
      UKBF2 =  -1.4546960D-5
      UKBF3 =   1.9047810D-6
      UKBF4 =  -1.1846940D-7
      UKBF5 =   3.5161530D-9
      UKBF6 =  -4.0853520D-11

      UKBDM3 =  0.0D0
      UKBDM2 =  1.0855000D-5
      UKBDM1 =  2.8391800D-4
      UKBD0 =  -0.048532D0
      UKBD1 =  -4.3289400D-5
      UKBD2 =   1.4902690D-5
      UKBD3 =  -1.9351350D-6
      UKBD4 =   1.1951780D-7
      UKBD5 =  -3.5235550D-9
      UKBD6 =   4.0708920D-11

      UKCFM3 =  0.0D0
      UKCFM2 =  1.9389200D-5
      UKCFM1 =  3.7263000D-4
      UKCF0 =   0.0485356D0
      UKCF1 =   1.5650290D-5
      UKCF2 =  -7.0134400D-6
      UKCF3 =   1.1275010D-6
      UKCF4 =  -8.2173200D-8
      UKCF5 =   2.7713320D-9
      UKCF6 =  -3.5625590D-11

      UKCDM3 =  0.0D0
      UKCDM2 = -1.8159700D-5
      UKCDM1 = -3.6821400D-4
      UKCD0 =  -0.0484683D0
      UKCD1 =  -1.4804670D-5
      UKCD2 =   6.8025790D-6
      UKCD3 =  -1.1017970D-6
      UKCD4 =   8.0248420D-8
      UKCD5 =  -2.6963500D-9

      UKCD6 =   3.4569950D-11

C  K1BF :
       AK1(1)= UKBFM3/(p*p*p)+ UKBFM2/(p*p)+ UKBFM1/p + UKBF0+UKBF1*p 
     >  +UKBF2*p*p+UKBF3*p*p*p+UKBF4*p*p*p*p+UKBF5*p*p*p*p*p 
     >  +UKBF6*p*p*p*p*p*p 
c           write(*,*) ' p = ',p
c           write(*,*) ' k1bf = ',ak1(1)
c           stop
C  K1CD :
       AK1(2)= UKCDM3/(p*p*p)+ UKCDM2/(p*p)+ UKCDM1/p + UKCD0+UKCD1*p 
     >  +UKCD2*p*p+UKCD3*p*p*p+UKCD4*p*p*p*p+UKCD5*p*p*p*p*p 
     >  +UKCD6*p*p*p*p*p*p 

C  K1AF :
       AK1(3)= UKAFM3/(p*p*p)+ UKAFM2/(p*p)+ UKAFM1/p + UKAF0+UKAF1*p 
     >  +UKAF2*p*p+UKAF3*p*p*p+UKAF4*p*p*p*p+UKAF5*p*p*p*p*p 
     >  +UKAF6*p*p*p*p*p*p 

C  K1BD :
       AK1(4)= UKBDM3/(p*p*p)+ UKBDM2/(p*p)+ UKBDM1/p + UKBD0+UKBD1*p 
     >  +UKBD2*p*p+UKBD3*p*p*p+UKBD4*p*p*p*p+UKBD5*p*p*p*p*p 
     >  +UKBD6*p*p*p*p*p*p 

C  K1CF :
       AK1(5)= UKCFM3/(p*p*p)+ UKCFM2/(p*p)+ UKCFM1/p + UKCF0+UKCF1*p 
     >  +UKCF2*p*p+UKCF3*p*p*p+UKCF4*p*p*p*p+UKCF5*p*p*p*p*p 
     >  +UKCF6*p*p*p*p*p*p 

C  K1AD :
       AK1(6)= UKADM3/(p*p*p)+ UKADM2/(p*p)+ UKADM1/p + UKAD0+UKAD1*p 
     >  +UKAD2*p*p+UKAD3*p*p*p+UKAD4*p*p*p*p+UKAD5*p*p*p*p*p  
     >  +UKAD6*p*p*p*p*p*p   

 

!!4.2 Main Magnet Sextupole strength vs momentum
!===========================

      TKADM3 = -4.54103D-5
      TKADM2 =  3.86648D-4
      TKADM1 =  -5.15221D-3
      TKAD0  = -6.23676D-3
      TKAD1  =  -8.21074D-5
      TKAD2  =  2.94841D-5
      TKAD3  =  -2.63597D-6  
      TKAD4  =  2.17817D-9
      TKAD5  =  6.02362D-9
      TKAD6  = -1.60702D-10

      TKAFM3 =  -2.11163D-5
      TKAFM2 =  2.31252D-4
      TKAFM1 =  -4.98909D-3
      TKAF0  =  -6.35613D-3
      TKAF1  =  -1.3545D-4
      TKAF2  =  5.20196D-5
      TKAF3  =  -5.93495D-6
      TKAF4  =  2.12422D-7
      TKAF5  =  6.36497D-11
      TKAF6  =  -9.88187D-11

      TKBDM3 =  -5.33308D-5  
      TKBDM2 =  4.61419D-4
      TKBDM1 =  -5.36794D-3
      TKBD0  =  -7.54801D-3
      TKBD1  =  -8.46257D-5
      TKBD2  =  3.25779D-5
      TKBD3  =  -3.46201D-6
      TKBD4  =  7.60854D-8
      TKBD5  =  3.2396D-9
      TKBD6  =  -1.24659D-10

      TKBFM3 =  -3.35980D-5
      TKBFM2 =  3.42004D-4
      TKBFM1 =  -5.26172D-3
      TKBF0  =  -7.65876D-3
      TKBF1  =  -1.20210D-4
      TKBF2  =  4.83680D-5
      TKBF3  =  -5.85176D-6
      TKBF4  =  2.31758D-7
      TKBF5  =  -1.23709D-9
      TKBF6  =  -7.77011D-11

      TKCDM3 =  5.46173D-5
      TKCDM2 =  -3.1184D-4
      TKCDM1 =  4.53022D-3
      TKCD0  =  -1.03323D-2
      TKCD1  =  -5.10398D-4
      TKCD2  =  1.74831D-4
      TKCD3  =  -1.90748D-5
      TKCD4  =  9.02456D-7
      TKCD5  =  -1.75112D-8
      TKCD6  =  7.24016D-11

      TKCFM3 =  3.24138D-5
      TKCFM2 =  -1.49126D-4
      TKCFM1 =  4.2301D-3
      TKCF0  =  -1.03626D-2
      TKCF1  =  -4.55906D-4
      TKCF2  =  1.55926D-4
      TKCF3  =  -1.70233D-5
      TKCF4  =  7.91965D-7
      TKCF5  =  -1.44494D-8
      TKCF6  =  3.8317D-11

!! Correction for B-dot effects is added in here.  K.Brown (8/31/98)
      TKDBDCOEF = -0.0025d0
      TKFBDCOEF = -0.0025d0

C  K2AD :
      AK2(6)= TKADM3/(p*p*p)+ TKADM2/(p*p)+ TKADM1/p + TKAD0+TKAD1*p 
     >  +TKAD2*p*p+TKAD3*p*p*p+TKAD4*p*p*p*p+TKAD5*p*p*p*p*p  
     >  +TKAD6*p*p*p*p*p*p - TKDBDCOEF*BDOT*ALA/ALA
      
C  K2AF :
      AK2(3)= TKAFM3/(p*p*p)+ TKAFM2/(p*p)+ TKAFM1/p + TKAF0+TKAF1*p 
     >  +TKAF2*p*p+TKAF3*p*p*p+TKAF4*p*p*p*p+TKAF5*p*p*p*p*p  
     >  +TKAF6*p*p*p*p*p*p - TKFBDCOEF*BDOT*ALA/ALA

C  K2BD :
      AK2(4)= TKBDM3/(p*p*p)+ TKBDM2/(p*p)+ TKBDM1/p + TKBD0+TKBD1*p 
     >  +TKBD2*p*p+TKBD3*p*p*p+TKBD4*p*p*p*p+TKBD5*p*p*p*p*p  
     >  +TKBD6*p*p*p*p*p*p - TKDBDCOEF*BDOT*ALB/ALA

C  K2BF :
      AK2(1)= TKBFM3/(p*p*p)+ TKBFM2/(p*p)+ TKBFM1/p + TKBF0+TKBF1*p 
     >  +TKBF2*p*p+TKBF3*p*p*p+TKBF4*p*p*p*p+TKBF5*p*p*p*p*p  
     >  +TKBF6*p*p*p*p*p*p - TKFBDCOEF*BDOT*ALB/ALA

C  K2CD :
      AK2(2)= TKCDM3/(p*p*p)+ TKCDM2/(p*p)+ TKCDM1/p + TKCD0+TKCD1*p 
     >  +TKCD2*p*p+TKCD3*p*p*p+TKCD4*p*p*p*p+TKCD5*p*p*p*p*p  
     >  +TKCD6*p*p*p*p*p*p - TKDBDCOEF*BDOT*ALC/ALA

C  K2CF :
      AK2(5)= TKCFM3/(p*p*p)+ TKCFM2/(p*p)+ TKCFM1/p + TKCF0+TKCF1*p 
     >  +TKCF2*p*p+TKCF3*p*p*p+TKCF4*p*p*p*p+TKCF5*p*p*p*p*p  
     >  +TKCF6*p*p*p*p*p*p - TKFBDCOEF*BDOT*ALC/ALA

c        write(89,*) ' p : ',p
c        write(89,*) ' madzg ak1 : ',ak1
c        write(89,*) ' madzg ak2 : ',ak2/2.d0
       return
       end
      FUNCTION STRCON(STR,STR2,
     >                         IS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL STRCON
      CHARACTER STR*(*), STR2*(*)
C     ---------------------------------------------------------------
C     .TRUE. if the string STR contains the string STR2 at least once
C     IS = position of first occurence of STR2 in STR 
C     (i.e.,STR(IS:IS+LEN(STR2)-1)=STR2)
C     ---------------------------------------------------------------
      INTEGER DEBSTR,FINSTR
      LNG2 = LEN(STR2(DEBSTR(STR2):FINSTR(STR2)))
      IF(LEN(STR).LT.LNG2) GOTO 1
      DO I = DEBSTR(STR), FINSTR(STR)-LNG2+1
        IF( STR(I:I+LNG2-1) .EQ. STR2 ) THEN
          IS = I 
          STRCON = .TRUE.
          RETURN
        ENDIF
      ENDDO
 1    CONTINUE
      STRCON = .FALSE.
      RETURN
      END
