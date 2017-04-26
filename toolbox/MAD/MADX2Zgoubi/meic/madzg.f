C23456789012345678901234567890123456789012345678901234567890123456789012
      implicit double precision (a-h,o-z)
      character fnr*80, fnw*80,  txt*132
      character txt590*590,keyword*20,txt20*20,txt60*60
      parameter(lr=9, lw=10, lout=14)
      parameter (mxl=15000)
      dimension iq(mxl)
      character kley*4, name*20
      parameter(mxk=20)
      character kle(mxk)*4, ny*1
      character(35) warn
      character(8) prtcl
      parameter (mss=3)
      character(20) stra(mss)
      logical TOMANY, okV

      logical ok, strcon
      integer debstr, finstr

      data kle/ 'DRIF', 'RBEN', 'SBEN', 'QUAD', 'SEXT', 'OCTU', 'MULT',
     >'SROT', 'YROT', 'MARK', 'KICK', 'HMON', 'VMON','HKIC', 'VKIC', 
     >'SOLE', 'RCOL', 'MATR', 'MONI', 'RFCA' /

C      write(*,*) ' name of the survey file to be translated:'
C      read(*,fmt='(A)') fnr
C      write(*,*) ' name of the zgoubi data file:'
C      read(*,fmt='(A)') fnw

      fnw='zgoubi.dat'

      fnr='madzg.in'
      if(fnr.ne.'madzg.in') then
        write(*,*)
        write(*,*) 'fnr MUST be  ''madzg.in'' '
        write(*,*)
      endif

      OPEN(UNIT=lout,FILE='madzg.out',err=97)
      OPEN(UNIT=lr,FILE=fnr,status='old',err=97)
 80   OPEN(UNIT=lw,FILE=fnw,STATUS='NEW',IOSTAT=IOS)
      IF(IOS .NE. 0) THEN
        CLOSE(lw,STATUS='DELETE')
        OPEN(UNIT=lw,FILE=fnw,STATUS='OLD')
        CLOSE(lw,STATUS='DELETE')
        GOTO 80
      ENDIF
      
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

c      write(*,*)
c      write(*,*) ' Convert DRIFT to zero-field MULTIPOLE (n/y) :'
c 6    read(*,*,err=6) ny
              ny = 'n'
      if(ny.ne.'y') ny='n' 
c      write(*,*) '               option set to ',ny
      call lmndri(ny)

c      write(*,*)
c      write(*,*) ' Translate SBEND to BEND or to MULTIPOLE (1/2) :'
c      read(*,*,err=7) kbm 
C            kbm = 1
            kbm = 1
C 7    continue
      if( kbm .ne. 2   ) kbm = 1
      call lmnkbm(kbm)
      if(kbm.eq.1   ) write(*,*) ' BEND will be used for SBEND'
      if(kbm.eq.2) write(*,*) ' MULTIPOLE will be used for SBEND'

C      okV = .false. 
      okV = .true. 
      write(*,*) ' Vertical bend :  ',okV
      call lmntV(okV)

c      write(*,*)
c      write(*,*) ' Translate RBEND to BEND or to MULTIPOLE (1/2) :'
c      read(*,*,err=71) kbmr 
            kbmr = 2
C 7    continue
      if( kbmr .ne. 2   ) kbmr = 1
      call lmnkbmr(kbmr)
      if(kbmr.eq.1   ) write(*,*) ' BEND will be used for RBEND'
      if(kbmr.eq.2) write(*,*) ' MULTIPOLE will be used for RBEND'

C Get momentum
 31   continue
        read(lr,fmt='(A)',end=85,err=99) txt590
        ok =       strcon(txt590,'@ PARTICLE',
     >                                        IS) 
      if(.not. ok) goto 31
      read(txt590(finstr(txt590)-9:finstr(txt590)),*) prtcl

 32   continue
        read(lr,fmt='(A)',end=85,err=99) txt590
        ok =       strcon(txt590,'@ PC',
     >                                  IS) 
      if(.not. ok) goto 32        
      read(txt590(finstr(txt590)-17:finstr(txt590)),*) pmom

c       write(*,*) pmom

      pmom =  1.d0 * 2.99792458d8 / 1.d9
      pmom =  60d9 / 1.d9
       
      it = 0
      call objet(lw,it,pmom*1.d9,prtcl,
     >                           bro)        
      bro = 1.d0 

        write(*,*) ' bro, p : ',bro, pmom

C Swallow rest of top of madzg.in file
 33   continue
        read(lr,fmt='(A)',end=85,err=99) txt590
c        write(*,*) txt590
        ok =       strcon(txt590,'KEYWORD',
     >                                    IS) 
      if(.not. ok) goto 33
      read(lr,fmt='(A)',end=85,err=99) txt590
c      write(*,*) txt590

      write(*,*) 
      write(*,*) ' Now translating. Busy...'

      ir = 0
      noel = 0
      sxl = 0.d0
      s1 = 0.d0
      s = 0.d0
 86   continue
      
        read(lr,fmt='(a)',end=85,err=99) txt590
        write(*,*) txt590

        read(txt590,fmt='(24e19.0,2e19.0,a)',end=85,err=99)
     >  xl, ang, ak1,ak2,ak3,ak4, tilt, e1,e2, h1,h2,sCntr,
     >  alfx,betx, alfy,bety, xco,yco, Dx,Dxp, Dy,Dyp, xmu,ymu, 
     >  hkic,vkic,
     >  txt60
C       txt60 includes : name, kley, solK

c        write(*,*)
c     >  xl, ang, ak1,ak2,ak3,ak4, tilt, e1,e2, h1,h2,sCntr,
c     >  alfx,betx, alfy,bety, xco,yco, Dx,Dxp, Dy,Dyp, xmu,ymu, txt60,
c     >  hkic,vkic

          
        call strget(txt60,mss,
     >                        nstr,stra)
        read(stra(1),*) name
        read(stra(2),*) kley
        solK = 0.d0
        if(nstr.eq.3) read(stra(3),*) solK

C        read(txt60,*,end=85,err=99) txt20,keyword
C        kley = keyword
C        name = txt20
c        write(*,*) '-----------------------------------'
c        write(*,*) kley, '  **', stra(1),'**',stra(2),'**',stra(3),'**'
c        write(*,*) kley, '  **', kley,'**',name,'**'

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
     >    (lw,lout,bro,frf,ik,noel,kley,name,xl,ang,ak1,ak2,ak3,ak4,
     >                                            tilt,e1,e2,h1,h2,
     >                                      hkic, vkic, solK,
     >                                      it)
c        call lmnt
c     >    (lw,bro,frf,ik,noel,kley,name,xl,ang,ak1,ak2,tilt,e1,e2,h1,h2,
c     >                                      it)
        ds = s - s1
c        if(abs(ds-xl).gt.1d-3) then
c          warn = 'WARNING : ds .ne. xl ; xl set to ds'
c          xl = ds
c          name = kley//' '//name
c          kley = 'DRIF'
c          it = it-1
c          call lmnt
c     >    (lw,bro,frf,ik,noel,kley,name,xl,ang,ak1,ak2,ak3,ak4,
c     >                                            tilt,e1,e2,h1,h2,
c     >                                      it)
cc          call lmnt
cc     >    (lw,bro,frf,ik,noel,kley,name,xl,ang,ak1,ak2,tilt,e1,e2,h1,h2,
cc     >                                      it)
c        else
c          warn = ' '
c        endif
        s1 = s
        sxl = sxl + xl
c        write(88,fmt='(1p,5(e16.8,1x),2(A16,1x),A35)') 
c     >                         s,sxl,s-sxl,xl,ds,kley,name,warn

      goto 86

 85   continue
c      write(lw,fmt='(A,T111,I6)') '''TWISS''',it
c      write(lw,fmt='(A)')   '2  1. 1. '
c      write(lw,fmt='(A,T111,I6)') '''MATRIX''',it
c      write(lw,fmt='(A)')   '1  11 '
c      it = it + 1
c      write(lw,fmt='(A,T111,I6)') '''END''',it
C      it = it + 1

      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''OPTIONS''',it
      write(lw,fmt='(a)') ' 1 1 '
      write(lw,fmt='(a)') ' WRITE ON'
      it = it + 1
      write(lw,fmt='(A,T111,I6)')  '''FAISCEAU''',it
      it = it + 1
      write(lw,fmt='(A,T111,I6)')  '''MARKER''     #End' ,it
      it = it + 1
      write(lw,fmt='(A,T111,I6)')  '''SPNPRT'''
      it = it + 1
c      write(lw,fmt='(A,T111,I6)')  '''REBELOTE''',it
c      write(lw,fmt='(A)')   '999999 0.2 99'
c      write(lw,fmt='(A,T111,I6)')  '''MATRIX''',it
c      write(lw,fmt='(A)')   ' 1  11  coupled'
      write(lw,fmt='(A,T111,I6)')  '''TWISS''',it
      write(lw,fmt='(A)')   ' 2  1.  1.  coupled   PRINT'
      it = it + 1
      write(lw,fmt='(A,T111,I6)')  '''END''',it

      write(*,*) ' Read ',ir,' elements from madzg.in file ',fnr
      write(*,*) ' Translated to',it,
     >' elements into the zgoubi data file ', fnw
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
        
      subroutine objet(lw,it,pmom,prtcl,
     >                            bro)
      implicit double precision (a-h,o-z)
      character(*) prtcl
      integer debstr, finstr

 1    continue

        bro = pmom /2.99792458d8    ! T.m 
        am = 938.27203d6
        gamma = sqrt(pmom*pmom + am*am)/am
      write(*,*) ' Problem rigidity (T.m), momentum (eV/c), G.gamma : ',
     >   bro,pmom,1.79284735d0*gamma
          
      write(lw,*) 'Generated by MADX -> Zgoubi translator'
      write(lw,fmt='(A,T111,I6)') '''OBJET''',it
      write(lw,FMT='(F19.6,7x,2(2x,A,F17.4))') bro*1.d3,
     >'reference rigidity (kG.cm) = ',pmom, 
     >',  G.gamma = ',1.79284735d0*gamma

c      write(lw,fmt='(A)') '2' 
c      write(lw,fmt='(A)') ' 1  1 ' 
c      write(lw,fmt='(A)') '.0 .2 .0 .2 0.  1.  ''o''' 
c      write(lw,fmt='(A)') ' 1' 
      write(lw,fmt='(A)') '5' 
      write(lw,fmt='(A)') ' .001 .01  .001 .01  .0 .0001 ' 
      write(lw,fmt='(A)') '.0 .0 .0 .0 .0 1. '

      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''PARTICUL''',it
      if    (prtcl(debstr(prtcl):finstr(prtcl)) .eq. 'ELECTRON') then
       write(lw,FMT='(A)') 
     >              '0.51099892 1.60217653e-19 1.15965218076e-3 0. 0.'
        write(lw,FMT='(A)') '''SRLOSS'''
        write(lw,FMT='(A)') '    0    ! .srLoss  '
        write(lw,FMT='(A)') '    BEND'
        write(lw,FMT='(A)') '    1  123456'
      else
       write(lw,FMT='(A)') '9.3827203E+02 1.602176487E-19 1.7928474 0 0' 
      endif

      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''SPNTRK''',it
      write(lw,fmt='(A)') ' 3'
c      write(lw,fmt='(A)') ' 4'
c      write(lw,fmt='(A)') '   1. 0. 0.'
      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''FAISCEAU''',it
      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''FAISTORE''',it
      write(lw,fmt='(A)') ' b_zgoubi.fai' 
      write(lw,fmt='(A)') ' 1'
      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''FAISCEAU''',it

c      it = it + 1
c      write(lw,fmt='(A,T111,I6)') '''OPTICS''',it
c      write(lw,FMT='(A)') 
c     > ' 2   Print out transport coeffs to zgoubi_MATRIX_out' 

      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''SCALING''' 
      write(lw,fmt='(A)') '1  3 '
      write(lw,fmt='(A)') 'BEND' 
      write(lw,fmt='(A)') '   -1   ' 
      write(lw,fmt='(1P,E16.8)') BRO
      write(lw,fmt='(A)') '1  '
      write(lw,fmt='(A)') 'MULTIPOL'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8)') BRO
      write(lw,fmt='(A)') '1      '
      write(lw,fmt='(A)') 'SOLENOID   OFF'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(A)') ' 40.5'
      write(lw,fmt='(A)') '1      '

      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''OPTIONS''',it
      write(lw,fmt='(a)') ' 1 1 '
      write(lw,fmt='(a)') ' WRITE ON'
      write(lw,fmt='(a)') ' '

      return
      end 

      subroutine lmnt
     >(lw,lout,bro,frf,ik,noel,kley,name,xl,ang,ak1,ak2,ak3,ak4,
     >                                        tilt,e1,e2,h1,h2,
     >                                      hkic, vkic,solK,
     >                                it)
c      subroutine lmnt
c     >(lw,bro,frf,ik,noel,kley,name,xl,ang,ak1,ak2,tilt,e1,e2,h1,h2,
c     >                                                       it)
      implicit double precision (a-h,o-z)
      character kley*(*)

      character(*) name
      character ny*1

      character*80 txfd, txfq,fqlhc , fdlhc ,fqrec , fdrec,fd,fq
      character*80 txffd,txffq,ffqlhc,ffdlhc,ffqrec,ffdrec,ffd,ffq
      character*80 ffdmu, fdmu, ffqmua,ffqmus,fqmua
      save txfd,txfq,txffd, txffq
      logical drimul
      save drimul

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
      save kpos, kbm, kbmr
      logical strcon

      logical okV, okVi, Vbnd
      save okV

      pi = 4.d0*atan(1.d0)
      zero=0.d0
      txtfrm = '(A,T12,A,T22,A,T111,i6)'

* meic. FM/VM Apr. 2014
      mstp = 1
c      if(strcon(name,'FF', 
c     >                    is)) then
c        write(*,*) ' mstp set to 2 ; name : ',name
c           read(*,*)
c         mstp = 2      
c      endif
 
C----- 1    2    3    4    5    6    7    8    9    10   11   12   13
C----- DRIF RBEN SBEN QUAD SEXT OCTU MULT SROT YROT MARK KICK HMON VMON
C----- 14   15   16   17   18   19   20
C----- HKIC VKIC SOLE RCOL MATR MONI RFCAVITY
       goto (
     > 1,   2 ,  2 ,  4,   5,   6,   7,   8,   9,   10,  11,  12,  12,
     > 11,  11,  13,  1,  18,   12,  1) ik

       stop ' Element was not found. '

 1    continue
C----- DRIF      
      if(xl .eq. 0.d0) then
        write(lw,fmt=txtfrm)'''MARKER''',kley,name,it
        goto 99
      else
        if(drimul) then
          dum=1.d-20
          write(lw,fmt=txtfrm) '''MULTIPOL''',kley,name,it
          write(lw,fmt='(I1,A)') i0,'  .Drift'
          write(lw,fmt='(G12.6,G8.2,G16.10,9G7.1)')
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
          write(lw,fmt=txtfrm) '''DRIFT''',kley,name,it
          write(lw,fmt='(F14.6)') xl*cm
        endif
      endif
      goto 99

 2    continue
       if(kbm.eq.1) then
C-------- Translate SBEND to BEND
         if(ik.eq.3) then
C---------- SBEN -> BEND
c           if(strcon(name,'BDOGLE',
c     >                            is) ) then
c              write(*,*) ' strcon(name, BDOGLE'
cc                  read(*,*)
c                goto 77
c            endif
           goto 3
         endif
       endif
       if(kbmr.eq.1) then
C-------- Translate RBEND to BEND
         if(ik.eq.2) then
C---------- RBEN -> BEND
           goto 21
         endif
       endif

 77    continue
C----- R- or SBEN -> MULTIPOL
            
c             mstp = 2   ! step size hlaved in RBEN

       if(xl .eq. 0.d0) then
         write(lw,fmt=txtfrm) '''MARKER''',kley,name,it
         goto 99
       endif

!      xlmag = xl
!      ro=xlmag/(2.d0*sin(ang/2.d0))
!      xlarc = ro * ang
      if(ik.eq.2) then
C case RBEN (madx twiss writes xl as arc length)
        xlarc = xl  !!!/(2.d0*sin(ang/2.d0)) * ang
      elseif(ik.eq.3) then
C case SBEN
        xlarc = xl
      else
        stop ' sbr lmnt, no such option ik '
      endif
      if(ang.ne.0.d0) then
        ro = xlarc /ang
        xlmag = 2.d0 * ro * sin (ang/2.d0)
        b1 = bro / ro                              ! dipole field (T)
        if(kpos.eq.3) then
          ale = -0.5d0 * ang
          xce=0.d0
          yce=0.d0
        elseif(kpos.eq.4) then
         ale = -0.5d0 * ang * 180/pi
         yce = cm*ro*(1-cos(ale*pi/180))
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
       write(lw,fmt=txtfrm) '''CHANGREF''',kley,name,it
       write(lw,fmt='(F13.9,3X,F13.9,3X,F13.9)') 0.d0,0.d0,ale 
       it = it+1
       write(lw,fmt=txtfrm) '''CHANGREF''',kley,name,it
       write(lw,fmt='(F13.9,3X,F13.9,3X,F13.9)') 0.d0,yce,0.d0
      endif

      write(lw,fmt=txtfrm) '''MULTIPOL''',kley,name,it
      write(lw,fmt='(I1,A)') i0,'  .Dip'

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
        off = ( si/2.d0/sqrk/(1.d0-co)- 1.d0/ak1)*ang
c        off = ( si/2.d0/sqrk/(1.d0-co)- 1.d0/xlarc/ak1)*ang
!        b1 = b1 + ak1 * off * Bro
!        b1 =  ak1 *Bro * si/2.d0/sqrk/(1.d0-co)*ang
!        write(*,fmt='(a,4f16.10)') name,off,ro,bro / ro, b1
      endif
! quadrupole field at .1 m (T) = G * a = k * Bro * a =  kl/xl * Bro * a
      b2 = ak1/xl *bro *.1d0
! sextupole field at .1 m (T)
      b3 = ak2/xl * bro *(x10/100.)**2 / 2.d0
C         write(*,*) ak2, bro, b3, x10

      if(abs(tilt).lt.1d-10) then 
        fac = 1.d0
      else
           if(strcon(name,'BDOGLE',
     >                            is) ) then
               fac = 1.d0
           else
              fac = 1.d-20
            endif
      endif
      write(lw,fmt='(F11.6,F7.2,3F13.8,7F4.1)')
     >xlmag*cm, x10, fac*b1*10.d0, fac*b2*10.d0, fac* b3*10.d0,
     > x0,x0,x0,x0,x0,x0,x0

      txt=txffd
      if(frf .eq. 0.) txt='0. 0. '//txt
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
C      write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(1p,4e16.8,A)') 
     >  tilt,tilt,tilt,tilt,'  0. 0. 0. 0. 0. 0.'
      if(xlmag*cm-int(xlmag*cm).gt.0.5) then 
            istepdip=( int(xlmag*cm)+1 ) 
      else
            istepdip=int(xlmag*cm)
      endif

            istepdip=   istepdip * mstp

      if(istepdip.lt.10) then
        write(lw,fmt='(2A)') ' #30|9|30    Dip',name
      elseif(istepdip.ge.10.and.istepdip.lt.100) then
        write(lw,fmt='(A,I2,A,2X,A)') ' #30|',istepdip,'|30    Dip',name
      elseif(istepdip.ge.100.and.istepdip.lt.1000) then
        write(lw,fmt='(A,I3,A,2X,A)') ' #30|',istepdip,'|30    Dip',name
      elseif(istepdip.ge.1000) then
        write(lw,fmt='(2A)') ' #30|9999|30    Dip',name
      endif
      if(kpos.ne.4) then
        if(abs(tilt).lt.1d-10) then
           write(6,fmt='(1P,I2,3(1x,e18.10))')  kpos,xce,yce,ale
           write(lw,fmt='(1P,I2,3(1x,e18.10))')  kpos,xce,yce,ale
        else
           write(6,fmt='(A)') ' 1. 0. 0. 0. '
           write(lw,fmt='(A)')  ' 1. 0. 0. 0. '
        endif
      else
       write(lw,*)  1,0.d0,0.d0,0.d0
      endif

      if(kpos.eq.4) then
       kley = ' '
       name = ' '
       it = it+1
       write(lw,fmt=txtfrm) '''CHANGREF''',kley,name,it
       write(lw,fmt='(F13.9,3X,F13.9,3X,F13.9)') 0.d0,-yce,0.d0
       it = it+1
       write(lw,fmt=txtfrm) '''CHANGREF''',kley,name,it
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

      Vbnd = okV .and. tilt.gt.1d-6
      if(Vbnd) then
        it = it+1
        write(lw,fmt=txtfrm) '''TRAROT''',name,'_TRAROT_I',it
C        write(lw,fmt=txtfrm) '''TRAROT''',name,kley,it
        write(lw,fmt='(a,1p,e14.6,a)') '0. 0. 0. ',tilt,' 0. 0.'
      endif

      write(lw,fmt=txtfrm) '''BEND''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .Bend'
      ro=xl/2./sin(ang/2.)
      b=ang * bro / xl *10.  !bro/ro*t2kg
C      xxl=ro*ang*cm
      ro=ro*cm
      write(lw,fmt='(F14.7,a,F15.8)') xl*cm,' 0. ',abs(b) 
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
      step = xl
      write(lw,fmt='(1P,E12.4,A)') dble(nint(step)),'  Bend'
C      write(lw,fmt='(F12.6,A4)') step * 10.,'  Bend'
C      write(lw,fmt='(A,2X,A)') ' #10|10|10   Bend',name
      if(kpos.eq.3) then
C        write(lw,fmt='(A)') '3 0. 0. 0.'
        write(lw,fmt='(A,1p,e16.8)') '3 0. 0. ',ale
      else
        write(lw,fmt='(1P,I1,2X,3G18.10)') kpos,xce,yce,ale
      endif

      if(Vbnd) then
        it = it+1
        write(lw,fmt=txtfrm) '''TRAROT''',name,'_TRAROT_O',it
C        write(lw,fmt=txtfrm) '''TRAROT''',name,kley,it
        write(lw,fmt='(a,1p,e14.6,a)') '0. 0. 0. ',-tilt,' 0. 0.'
      endif
      goto 99

c 21   continue
cC----- RBEN -> BEND
c      if(xl * ang .eq. 0.d0) then
c        if(xl .eq. 0.d0) then
c          write(lw,fmt=txtfrm) '''MARKER''',kley,name,it
c        else
c          write(lw,fmt=txtfrm) '''DRIFT''',kley,name,it
c          write(lw,fmt='(1P,E14.6)') xl*cm
c        endif
c        goto 99
c      endif
c
c      Vbnd = okV .and. tilt.gt.1d-6
c      if(Vbnd) then
c        it = it+1
c        write(lw,fmt=txtfrm) '''TRAROT''',name,'_TRAROT_I',it
cC        write(lw,fmt=txtfrm) '''TRAROT''',name,kley,it
c        write(lw,fmt='(a,1p,e14.6,a)') '0. 0. 0. ',tilt,' 0. 0.'
c      endif
c
c      write(lw,fmt=txtfrm) '''BEND''',kley,name,it
c      write(lw,fmt='(I1,A)') i0,'  .Bend'
c      ro=xl/2./sin(ang/2.)
c      b=ang * bro / xl *10.  !bro/ro*t2kg
cC      xxl=ro*ang*cm
c      ro=ro*cm
c
cC err. FM June 2009      write(lw,fmt='(2F14.7,F15.8)') xl,tilt,b 
c      if(abs(tilt).lt.1d-10) then 
c        fac = 1.d0
c      else
c        fac = 1.d-20
c      endif
c      write(lw,fmt='(2F16.9,F15.8)') xl*cm,tilt,fac*b 
c      te=ang/2.+e1
c      if(frf .eq. 0.) then
c        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,te
c      else 
c        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,te
c      endif
c      write(lw,fmt='(A)') txfd
c      ts=ang/2.+e2
c      if(frf .eq. 0.) then
c        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,ts
c      else
c        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,ts
c      endif
c      write(lw,fmt='(A)') txfd
c      step = xl
c      write(lw,fmt='(1P,E12.4,A)') dble(nint(step)),'  Bend'
cC      write(lw,fmt='(F12.6,A4)') step * 10.,'  Bend'
cC      write(lw,fmt='(A,2X,A)') ' #10|10|10   Bend',name
c      if(kpos.eq.3) then
c        if(abs(tilt).lt.1d-10) then 
c          write(lw,fmt='(A)') '3 0. 0. 0.'
c        else
c          write(lw,fmt='(A)') '1 0. 0. 0.'
c        endif
c      else
c        write(lw,fmt='(1P,I2,3(1x,e18.10))') kpos,xce,yce,ale
c      endif
c      goto 99

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

      Vbnd = okV .and. tilt.gt.1d-6
      if(Vbnd) then
        it = it+1
        write(lw,fmt=txtfrm) '''TRAROT''',name,'_TRAROT_I',it
C        write(lw,fmt=txtfrm) '''TRAROT''',name,'_TRAROT_I',it
        write(lw,fmt='(a,1p,e14.6,a)') '0. 0. 0. ',tilt,' 0. 0.'
      endif

      ro=xl/ang
      b=bro/ro*t2kg

      if(b<0.d0) then
        write(lw,fmt=txtfrm) '''YMY''',name,'_YMY_I',it
        write(lw,fmt=txtfrm) '''BEND''',name,'NEG_B',it
      else
        write(lw,fmt=txtfrm) '''BEND''',name,kley,it
      endif
      write(lw,fmt='(I1,A)') i0,'  .Bend'
      xxl = 2.d0*ro*sin(xl/2.d0/ro)
      write(lw,fmt='(F14.7,a,F15.8)') xxl*cm,' 0. ',abs(b) 
      if(b<0.d0) then 
        te = -e1
        ale = ang/2.d0
      else
        te=e1
        ale = -ang/2.d0
      endif
      if(frf .eq. 0.) then
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,te
      else
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,te
      endif
      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
      ts=e2
      if(b<0.d0) ts = -ts
      if(frf .eq. 0.) then
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,ts
      else
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,ts
      endif
      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
C      istepdip=int(xxl*cm/3.0d0)
      istepdip=int(xxl*cm)
      if(istepdip .ge. 100) then
        write(lw,fmt='(A,I3,A,2X,A)') '#200|',istepdip,
     + '|200    Bend',name
      elseif(istepdip .ge. 10  .and. istepdip .lt. 100) then
        write(lw,fmt='(A,I2,A,2X,A)') '#200|',istepdip,
     + '|200    Bend',name
      elseif(istepdip .lt. 10) then
        write(lw,fmt='(2A)') '#200|5|200    Bend',name 
      endif
      if(kpos.eq.3) then
C        write(lw,fmt='(A)') '3 0. 0. 0.'
        write(lw,fmt='(A,1p,e16.8)') '3 0. 0. ',ale
      else
        write(lw,fmt='(1P,I1,2X,3G18.10)') kpos,xce,yce,ale
      endif
      if(b<0.d0) 
     >write(lw,fmt=txtfrm) '''YMY''',name,'_YMY_O',it

      if(Vbnd) then
        it = it+1
        write(lw,fmt=txtfrm) '''TRAROT''',name,'_TRAROT_O',it
C        write(lw,fmt=txtfrm) '''TRAROT''',name,kley,it
        write(lw,fmt='(a,1p,e14.6,a)') '0. 0. 0. ',-tilt,' 0. 0.'
      endif

      goto 99

c 3    continue
cC----- SBEN -> BEND
c      if(xl * ang .eq. 0.d0) then
c        if(xl .eq. 0.d0) then
c          write(lw,fmt=txtfrm) '''MARKER''',kley,name,it
c        else
c          write(lw,fmt=txtfrm) '''DRIFT''',kley,name,it
c          write(lw,fmt='(1P,E14.6)') xl*cm
c        endif
c        goto 99
c      endif
c
c      Vbnd = okV .and. tilt.gt.1d-6
c      if(Vbnd) then
c        it = it+1
c        write(lw,fmt=txtfrm) '''TRAROT''',name,'_TRAROT_I',it
cC        write(lw,fmt=txtfrm) '''TRAROT''',name,'_TRAROT_I',it
c        write(lw,fmt='(a,1p,e14.6,a)') '0. 0. 0. ',tilt,' 0. 0.'
c      endif
c
c      write(lw,fmt=txtfrm) '''BEND''',kley,name,it
c      write(lw,fmt='(I1,A)') i0,'  .Bend'
c      ro=xl/ang
c      xxl = 2.d0*ro*sin(ang/2.d0)
c      b=bro/ro*t2kg
c              write(lout,*) ' ********* tilt xl ',tilt,xxl
c      if(abs(tilt).lt.1d-10) then 
c        fac = 1.d0
c      else
c        fac = 1.d-20
c      endif
cc      write(lw,fmt='(1P,3(1X,E17.9))') xxl*cm,tilt,fac*b 
c      write(lw,fmt='(1P,3(1X,E17.9))') xxl*cm,tilt,b 
c      te=e1
c      if(frf .eq. 0.) then
c        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,te
c      else
c        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,te
c      endif
c      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
c      ts=e2
c      if(frf .eq. 0.) then
c        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,ts
c      else
c        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,ts
c      endif
c      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
cC      istepdip=int(xxl*cm/3.0d0)
c      istepdip=int(xxl*cm) 
c
C             mstp = 2   ! step size hlaved in SBEN
c            istepdip=   istepdip * mstp
c
c      if(istepdip .ge. 100) then
c        write(lw,fmt='(A,I3,A,2X,A)') '#200|',istepdip,
c     + '|200    Bend',name
c      elseif(istepdip .ge. 10  .and. istepdip .lt. 100) then
c        write(lw,fmt='(A,I2,A,2X,A)') '#200|',istepdip,
c     + '|200    Bend',name
c      elseif(istepdip .lt. 10) then
c        write(lw,fmt='(2A)') '#200|5|200    Bend',name 
c      endif
c      if(kpos.eq.3) then
c        if(abs(tilt).lt.1d-10) then 
c          write(lw,fmt='(A,1p,2x,e14.6)') '3 0. 0. ',-ang/2.d0
c        else
c          write(lw,fmt='(A)') '1 0. 0. 0.'
c        endif
c      else
c        write(lw,fmt='(1P,I2,3(1x,e18.10))') kpos,xce,yce,ale
c      endif
c      goto 99

 4    continue
C----- QUAD
      if(xl .eq. 0.d0) then
        write(lw,fmt=txtfrm) '''MARKER''',kley,name,it
        goto 99
      endif
      write(lw,fmt=txtfrm) '''MULTIPOL''',kley,name,it
      write(lw,fmt='(I1,A)') i0,'  .Quad'
      write(lw,fmt='(F12.6,F6.2,2F16.10,8F4.1)')
     >xl*cm,x10,x0,ak1/xl*bro,x0,x0,x0,x0,x0,x0,x0,x0
      txt = txffq
      txt=txt//name
      if(frf .eq. 0.) txt = '0. 0. '//txt
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfq
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfq
      write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      if(xl*cm-int(xl*cm).gt.0.5) then 
            istepdip=int(xl*cm)+1
      else
            istepdip=int(xl*cm)
      endif

             mstp = 4   ! step size hlaved in QUAD
      istepdip = istepdip * mstp

      if(istepdip.lt.6) then
        write(lw,fmt='(2A)') ' #30|5|30   Quad',name
      elseif(istepdip.ge.10.and.istepdip.lt.100) then
        write(lw,fmt='(A,I2,A,2X,A)') ' #30|',istepdip,'|30   Quad',name
      elseif(istepdip.ge.100.and.istepdip.lt.1000) then
        write(lw,fmt='(A,I3,A,2X,A)') ' #30|',istepdip,'|30   Quad',name
      elseif(istepdip.ge.1000) then
        write(lw,fmt='(2A)') ' #30|9999|30   Quad',name
      endif
c      write(lw,fmt='(A,2X,A)') '#10|10|10  Quad',name
      write(lw,fmt='(A)') '1 0. 0. 0.'
      goto 99

 5    continue
C----- SEXT ***
      if(xl .eq. 0.d0) then
        write(lw,fmt=txtfrm) '''MARKER''',kley,name,it
        goto 99
      endif
      write(lw,fmt=txtfrm) '''MULTIPOL''',kley,name,it
      write(lw,fmt='(I1,A)') i0,'  .Sext'
      b3 = 10.d0 * ak2/xl * bro *(x10/100.)**2 / 2.d0 
      write(lout,*) ' SEXT         b3 = ',b3,' kG,  xl = ',xl
      write(lw,fmt='(F12.6,F6.2,2F10.4,F16.10,7F4.1)')
     >xl*cm,x10,x0,x0,b3,x0,x0,x0,x0,x0,x0,x0
      write(lw,fmt='(A)') '0.00 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') '0.00 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      if(xl*cm-int(xl*cm).gt.0.5) then 
            istepdip=int(xl*cm)+1
      else
            istepdip=int(xl*cm)
      endif
      if(istepdip.lt.10) then
        write(lw,fmt='(2A)') ' #30|9|30   Sext',name
      elseif(istepdip.ge.10.and.istepdip.lt.100) then
        write(lw,fmt='(A,I2,A,2X,A)') ' #30|',istepdip,'|30   Sext',name
      elseif(istepdip.ge.100.and.istepdip.lt.1000) then
        write(lw,fmt='(A,I3,A,2X,A)') ' #30|',istepdip,'|30   Sext',name
      elseif(istepdip.ge.1000) then
        write(lw,fmt='(2A)') ' #30|9999|30   Sext',name
      endif
C      write(lw,fmt='(A,2X,A)') '#10|10|10  Sext',name
      write(lw,fmt='(A)') '1 0. 0. 0.'
      goto 99

 6    continue
C----- OCTU ***
      if(xl .eq. 0.d0) then
        write(lw,fmt=txtfrm) '''MARKER''',kley,name,it
C        it = it - 1
        goto 99
      endif
      write(lw,fmt=txtfrm) '''MULTIPOL''',kley,name,it
      write(lw,fmt='(I1,A)') i0,'  .Octu'
      write(lw,fmt='(F12.6,F6.2,2F16.10,4F9.4)')
     >xl*cm,x10,x0,x0,x0,x0,x0,x0
      write(lw,fmt='(A)') '0.000 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') '0.000 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      if(xl*cm-int(xl*cm).gt.0.5) then 
            istepdip=int(xl*cm)+1
      else
            istepdip=int(xl*cm)
      endif
      if(istepdip.lt.10) then
        write(lw,fmt='(2A)') ' #30|9|30   Octu',name
      elseif(istepdip.ge.10.and.istepdip.lt.100) then
        write(lw,fmt='(A,I2,A,2X,A)') ' #30|',istepdip,'|30   Octu',name
      elseif(istepdip.ge.100.and.istepdip.lt.1000) then
        write(lw,fmt='(A,I3,A,2X,A)') ' #30|',istepdip,'|30   Octu',name
      elseif(istepdip.ge.1000) then
        write(lw,fmt='(2A)') ' #30|9999|30   Octu',name
      endif
C      write(lw,fmt='(A,2X,A)') '#20|20|20  Octu',name
      write(lw,fmt='(A)') '1 0. 0. 0.'
       goto 99

 7    continue
C----- MULT
      if(abs(ak1)+abs(ak2)+abs(ak3)+abs(ak4) .gt. 1.d-30) then
        xl = 1.d-5
      else
        write(lw,fmt=txtfrm)'''MARKER''',kley,name,it
        goto 99
      endif
! quadrupole field at .1 m (kG)
      b2 = ak1/xl * bro   
! sextupole field at .1 m (kG)
      b3 = 10.d0 * ak2/xl * bro *(x10/100.)**2 / 2.d0  
! octupole field at .1 m (kG)
      b4 = 0.d0
      b5 = 0.d0
      b6 = 0.d0
      write(lw,fmt=txtfrm) '''MULTIPOL''',kley,name,it
      write(lw,fmt='(I1,A)') i0,'  .Mult'
      write(lw,fmt='(F12.6,F6.2,2F16.10,8F4.1)')
     >xl*cm,x10,x0,b2,b3,b4,b5,b6,x0,x0,x0,x0
c      write(*,*)
c     >b1,b2,b3,b4,b5,b6,bro,' b1,b2,b3,b4,b5,b6, bro'
c          stop
c      write(*,*) ak1,' ak1 '
c      write(87,*)
c     > abs(ak1)+abs(ak2)+abs(ak3)+abs(ak4),' abs(k1)+abs(k2...'
      write(lw,fmt='(A)') '00.00 00.0  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') '00.00 00.0  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      if(xl*cm-int(xl*cm).gt.0.5) then 
            istepdip=int(xl*cm)+1
      else
            istepdip=int(xl*cm)
      endif

      if(istepdip.lt.10) then
        write(lw,fmt='(2A)') ' #30|9|30   Mult',name
      elseif(istepdip.ge.10.and.istepdip.lt.100) then
        write(lw,fmt='(A,I2,A,2X,A)') ' #30|',istepdip,'|30   Mult',name
      elseif(istepdip.ge.100.and.istepdip.lt.1000) then
        write(lw,fmt='(A,I3,A,2X,A)') ' #30|',istepdip,'|30   Mult',name
      elseif(istepdip.ge.1000) then
        write(lw,fmt='(2A)') ' #30|9999|30   Mult',name
      endif
C      write(lw,fmt='(A,2X,A)') '#20|20|20  Mult',name
      write(lw,fmt='(A)') '1 0. 0. 0.'
       goto 99

 8    continue
C----- SROT ***
      write(lw,fmt=txtfrm) '''TRAROT''',kley,name,it
      goto 99

 9    continue
C----- YROT ***
      write(lw,fmt=txtfrm) '''TRAROT''',kley,name,it
      goto 99

 10   continue
C----- MARK ***
      write(lw,fmt=txtfrm) '''MARKER''',kley,name,it
C      it = it - 1
      goto 99

 11   continue
C----- KICK ***
      xxl = xl
      if(xl .eq. 0.d0) xxl = 1.e-6
      if(ik.eq.14) then
c--------- hkicker
        diptlt = 0.d0
        b = -bro*hkic/xxl*t2kg
      elseif(ik.eq.15) then
c--------- vkicker
        diptlt = pi/2.d0
        b = -bro*vkic/xxl*t2kg
      endif
c             write(*,*) ' KICK, hkic, vkic :', hkic, vkic
      write(lw,fmt=txtfrm) '''MULTIPOL''',kley,name,it
      write(lw,fmt='(I1,A)') i0,'  .kicker'
      write(lw,fmt='(E14.6,F6.2,E14.6,F12.6,8F4.1)')
     >xxl*cm,x10,b,x0,x0,x0,x0,x0,x0,x0,x0,x0
      txt = '.0 .0  1.00 0.00 0.00 0.00 0.00 0. 0. 0. 0.'
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(f12.9,A)') diptlt,' 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(A)') '#20|4|20  Kick'
      write(lw,fmt='(A)') '1 0. 0. 0.'
      goto 99

 12   continue
C----- HMON, VMON
      write(lw,fmt=txtfrm) '''DRIFT''',kley,name,it
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
        write(lw,fmt=txtfrm) '''MARKER''',kley,name,it
C        it = it - 1
        goto 99
      endif
C      r0 = 111.11
      r0 = 2.
      ffxtent= 2.5 * r0
      write(lw,fmt=txtfrm) '''SOLENOID''','  OFF ',kley,it
C      write(lw,fmt=txtfrm) '''SOLENOID''',kley,name,it
      write(lw,fmt='(I1,A)') i0,'  .soleno'
      write(lw,fmt='(1P,G14.7,G10.4,G14.7,a)') 
     >                     xl*cm,r0,solK*bro*t2kg,' 1'
      write(lw,fmt='(2F10.2)') ffxtent, ffxtent
C      write(lw,fmt='(A,2X,A)') '1.  Soleno',name
      write(lw,fmt='(A,2X,A)') '1.  Soleno',name
      write(lw,fmt='(A)') '1 0. 0. 0. '
      goto 99

 18   continue
C----- MATR
      if(xl .eq. 0.d0) then
        write(lw,fmt=txtfrm)'''MARKER''',kley,name,it
        goto 99
      else
        stop ' Element MATR not translated' 
      endif
      goto 99

 99    continue
       write(6,fmt='(A,T12,A)') kley,name
       return

      entry lmntV(okVi)
      okV = okVi
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

      entry lmnkbm(kbmi)
      kbm =kbmi 
      return

      entry lmnkbmr(kbmri)
      kbmr =kbmri 
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
      SUBROUTINE STRGET(STR,MSS,
     >                          NST,STRA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*) STRA(*)
      CHARACTER(*) STR
C     ------------------------------------------------------
C     Extract substrings #1 up to #MSS, out of string STR. 
C     Strings are assumed spaced by (at least) one blank. 
C     They are saved in  array STRA, and their total number 
C     (possibly < mss) is NST.
C     ------------------------------------------------------
      INTEGER DEBSTR, FINSTR

      CHARACTER(2000) STR0

      IF(LEN(STR0).LT.LEN(STR)) 
     >stop ' Sbr strget : increase length of string str0.'

      STR0 = STR
      IE = FINSTR(STR)
      NST = 0
      I2 = 1

 1    CONTINUE

        IF(STR(I2:I2) .EQ. ' '  .OR. 
     >     STR(I2:I2) .EQ. ',') THEN
          I2 = I2 + 1
          IF(I2 .LE. IE) GOTO 1
        ELSE
          I1 = I2
 2        CONTINUE
          I2 = I2 + 1
          IF(I2 .LE. IE) THEN
            IF(STR(I2:I2) .EQ. ' '  .OR. 
     >         STR(I2:I2) .EQ. ',') THEN
              IF(NST .LT. MSS) THEN
                NST = NST + 1
                STRA(NST) = STR(I1:I2-1)
                I2 = I2 + 1
                GOTO 1
              ENDIF
            ELSE
              GOTO 2
            ENDIF
          ELSE
            IF(STR(I2-1:I2-1) .NE. ' ' .AND.
     >         STR(I2-1:I2-1) .NE. ',') THEN
              IF(NST .LT. MSS) THEN
                NST = NST + 1
                STRA(NST) = STR(I1:I2-1)
              ENDIF
            ENDIF
          ENDIF
        ENDIF

      STR = STR0

c      call ZGNOEL(
c     >             NOEL)
c       if(noel.eq.89) then
c           write(*,*) ' strget  //////////////////'
c           write(*,*) ' strget  NST = ', nst
c           write(*,*) ' strget ', (stra(i),i=1,nst)
c            write(*,*) ' strget  //////////////////'
c              read(*,*)     
c          endif

      RETURN
      END
