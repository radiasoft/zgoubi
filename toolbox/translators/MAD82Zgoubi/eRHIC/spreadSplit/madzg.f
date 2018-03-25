C23456789012345678901234567890123456789012345678901234567890123456789012
      implicit double precision (a-h,o-z)
      character fnr*80, fnw*80, fnw2*81, txt*132, titl*132
      parameter(lr=9, lw=10, lw2=12, lout=14)
      parameter (mxl=15000)
      dimension iq(mxl)
      character kley*4, name*16
      parameter(mxk=20)
      character kle(mxk)*4, ny*1
      character*35 warn
      logical TOMANY, okV

      data kle/ 'DRIF', 'RBEN', 'SBEN', 'QUAD', 'SEXT', 'OCTU', 'MULT',
     >'SROT', 'YROT', 'MARK', 'KICK', 'HMON', 'VMON','HKIC', 'VKIC', 
     >'SOLE', 'RCOL', 'MATR', 'INST' , 'MONI' /

C      write(*,*) ' name of the survey file to be translated:'
C      read(*,fmt='(A)') fnr
C      write(*,*) ' name of the zgoubi data file:'
C      read(*,fmt='(A)') fnw

      fnr='twiss'
c      fnr='survey'
      if(fnr.ne.'twiss' .and. fnr.ne.'survey') then
        write(*,*)
        write(*,*) 'fnr MUST be either ''twiss'' or ''survey'''
        write(*,*)
      endif
      fnw='trad.out'
      fnw2='trad.out2'

      OPEN(UNIT=lout,FILE='madzg.out',err=97)
      OPEN(UNIT=lr,FILE=fnr,status='old',err=97)
 80   OPEN(UNIT=lw,FILE=fnw,STATUS='NEW',IOSTAT=IOS)
      IF(IOS .NE. 0) THEN
        CLOSE(lw,STATUS='DELETE')
        OPEN(UNIT=lw,FILE=fnw,STATUS='OLD')
        CLOSE(lw,STATUS='DELETE')
        GOTO 80
      ENDIF
 82   OPEN(UNIT=lw2,FILE=fnw2,STATUS='NEW',IOSTAT=IOS)
      IF(IOS .NE. 0) THEN
        CLOSE(lw2,STATUS='DELETE')
        OPEN(UNIT=lw2,FILE=fnw2,STATUS='OLD')
        CLOSE(lw2,STATUS='DELETE')
        GOTO 82
      ENDIF
      
      read(lr,fmt='(A)',err=99) txt    
      read(lr,fmt='(A)',err=99) titl
      write(*,*) txt
      write(lw,fmt='(A60)') titl

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
      write(*,*) ' Convert DRIFT to zero-field MULTIPOLE (n/y) :'
c 6    read(*,*,err=6) ny
              ny = 'n'
      if(ny.ne.'y') ny='n' 
      write(*,*) '               option set to ',ny
      call lmndri(ny)

      write(*,*)
      write(*,*) ' Translate S/RBEND to BEND or to MULTIPOLE (1/2) :'
C      read(*,*,err=7) kbm
            kbm = 1
 7    continue
      if( kbm .ne. 2   ) kbm = 1
      if(kbm.eq.1   ) write(*,*) ' BEND will be used for S/RBEND'
      if(kbm.eq.2) write(*,*) ' MULTIPOLE will be used for S/RBEND'
      call lmnkbm(kbm)

C      okV = .false. 
      okV = .true. 
      write(*,*) ' Vertical bend :  ',okV
      call lmntV(okV)

      if(fnr.eq.'twiss') i2 = 0
      if(fnr.eq.'survey') i2 = 4
      do i=1,i2
        read(lr,fmt='(A)',err=99) txt
      enddo        

      write(*,*) 
      write(*,*) ' Now translating. Busy...'

      if    (fnr.eq.'survey') then
        read(lr,fmt='(A4,A16,F12.6,3E16.9,/,5E16.9)',end=85,err=99)
     >   kley, name, xl, ang, ak1, ak2, tilt, e1, e2, h1, h2
        read(lr,fmt='(1P,4E16.9)',end=85,err=99) r1,r2,r3,s
        read(lr,fmt='(A)',end=85,err=99) txt
      elseif(fnr.eq.'twiss') then
        read(lr,fmt='(A4,A16,F12.6,3E16.9,/,5E16.9)',end=85,err=99)
     >   kley, name, xl, ang, ak1, ak2, tilt, e1, e2, h1, h2
C                                                
        read(lr,fmt='(1P,5E16.9)',end=85,err=99) alfx,betx,gamx,dx,dxp 
        read(lr,fmt='(1P,5E16.9)',end=85,err=99) alfy,bety,gamy,dy,dyp 
C                                                
        read(lr,fmt='(1P,5E16.9)',end=85,err=99) r21,r22,r23,r24,r25
        write(lw2,*) s,xl, ang, ak1, ak2, tilt, e1, e2, h1, h2
        write(lw2,*)  alfx,betx,gamx,dx,dxp ,' alfx,betx,gamx,dx,dxp '
        write(lw2,*)  alfy,bety,gamy,dy,dyp ,' alfy,bety,gamy,dy,dyp '
      endif

      it = 1
        write(*,*)  alfx,betx,gamx,dx,dxp ,' alfx,betx,gamx,dx,dxp '
        write(*,*)  alfy,bety,gamy,dy,dyp ,' alfy,bety,gamy,dy,dyp '
c              read(*,*)
      call objet(lw,it,bro, titl,
     > alfx,betx,dx,dxp, 
     > alfy,bety,dy,dyp )

      ir = 0
      noel = 0
      sxl = 0.d0
      s1 = 0.d0
      s = 0.d0
 86   continue
      
      if    (fnr.eq.'survey') then
        read(lr,fmt='(A4,A16,F12.6,3E16.9,/,5E16.9)',end=85,err=99)
     >   kley, name, xl, ang, ak1, ak2, tilt, e1, e2, h1, h2
        read(lr,fmt='(1P,4E16.9)',end=85,err=99) r1,r2,r3,s
        read(lr,fmt='(A)',end=85,err=99) txt
      elseif(fnr.eq.'twiss') then
        read(lr,fmt='(A4,A16,F12.6,3E16.9,/,5E16.9)',end=85,err=99)
     >   kley, name, xl, ang, ak1, ak2, tilt, e1, e2, h1, h2
C                                                
        read(lr,fmt='(1P,5E16.9)',end=85,err=99) alfx,betx,gamx,dx,dxp 
        read(lr,fmt='(1P,5E16.9)',end=85,err=99) alfy,bety,gamy,dy,dyp 
C                                                
        read(lr,fmt='(1P,5E16.9)',end=85,err=99) r21,r22,r23,r24,r25
        write(lw2,*) s,xl, ang, ak1, ak2, tilt, e1, e2, h1, h2
        write(lw2,*)  alfx,betx,gamx,dx,dxp ,' alfx,betx,gamx,dx,dxp '
        write(lw2,*)  alfy,bety,gamy,dy,dyp ,' alfy,bety,gamy,dy,dyp '
      endif

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
        write(88,fmt='(1p,5(e16.8,1x),2(A16,1x),A35)') 
     >                         s,sxl,s-sxl,xl,ds,kley,name,warn

      goto 86

 85   continue

      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''FAISCEAU''',it

      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''MATRIX''',it
      write(lw,fmt='(A)')   '1  0 '

      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''SPNPRT''',it
      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''SRPRNT''',it
c      it = it + 1
c      write(lw,fmt='(A,T111,I6)') '''TWISS''',it
c      write(lw,fmt='(A)')   '1  1. 1. '

c      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''END''',it

C      it = it + 1
      write(lw,fmt='(///,A)') ' '
      write(lw,fmt='(///,A)') '''REBELOTE'''
      write(lw,fmt='(A)')   '600 0.1 99'
      write(*,*) ' Read ',ir,' elements from the survey file ',fnr
      write(*,*) ' Translated to',it,
     >' elements into the zgoubi data file "trad.out"', fnw
      write(*,*) ' end of survey file '
      goto 999

 97   write(*,*) ' error open survey file'
      goto 999
 99   write(*,*) ' error during read of survey file'
      goto 999

 999  continue
      close(lr)
      close(lw)
      close(lw2)
      stop
      end
        
      subroutine objet(lw,it,bro, titl,
     > alfx,betx,dx,dxp, 
     > alfy,bety,dy,dyp )
      implicit double precision (a-h,o-z)
C      parameter(G = 1.79284735, am = 938.27203d6, c=2.99792458d8)
      character(*) titl
      parameter(G=1.15965218076e-3, am=0.51099892d6, c=2.99792458d8)
!      parameter(GgDef = 37.4448168746)   !  16.5GeV 
!      parameter(GgDef = 3.00011777488)   !  1.322GeV 
      logical strcon, ok, ok2
      character(20) stra(20)
      integer debstr, finstr

      ok=strcon(titl,'GeV',
     >                     is)
      if(ok) then
        call strget(titl,20,
     >                      NSS,stra)
         
        if(nss.gt.0) then
          ok2 = .false.
          i = 1
          dowhile ((i .le. nss) .and.  (ok2 .eqv. .false.))
            stra(i) = stra(i)(debstr(stra(i)):finstr(stra(i)))
            if(strcon(stra(i),'GeV',
     >                            is)) then
              ok2 = .true.
              if(is.gt.1) then
                read(stra(i)(1:is-1),*) TGeV
              else
                read(stra(i-1),*) TGeV
              endif
            endif
            i = i + 1
          enddo
        else
          TGeV = 1d-6
        endif
      else
        TGeV = 1d-6
      endif

      amGeV = am*1d-9
C      TGeV = 1.665d0
      pGeV = sqrt( (TGeV + amGeV)**2 - amGeV**2)
      bta = pGeV/sqrt(pGeV*pGeV + amGeV**2)
      gam = pGeV/bta / amGeV
      GgDef = G * gam
C      TGeV = (gam-1.d0)*(am*1d-9)

c 1    continue
C      write(*,*) 'Give  reference  G*gamma  (default is ) ',GgDef
        TDef = GgDef / G * am - am
      write(*,*) 'Give  reference kin-energy/eV (default is ) ',TDef
C      read(*,*,err=2,end=2) Gg
CCCCCCCCCCCCCCCC      read(*,*,err=2,end=2) Tr
CCCCCCCCCCCCCCC      goto 3
 2    continue
      Tr = TDef
C      Gg = GgDef
 3    continue
      
      broDef = sqrt(Tr*(Tr+2.d0*am))/c    ! T.m 

      write(*,*)  ' Reference rigidity :'
      write(*,FMT='(F19.6,A,/,A)') 
     >      broDef  !!!!!!,' T.m','  Enter desired value'
c      read(*,*,err=1) broDef

      write(*,*) 'Give working kin-energy (GeV) ',TGeV
C      read(*,*,err=2,end=2) Gg
CCCCCCCCCCCCCCCCC      read(*,*,err=2,end=2) TGeV
      T = TGeV*1d9 
      broW = sqrt(T*(T+2.d0*am))/c    ! T.m 

      bro = broW/broDef    
           write(*,*) 'bro = ',bro
c                 read(*,*)
               
      write(lw,fmt='(A,T111,I6)') '''SYSTEM'''
      write(lw,fmt='(A,T111,I6)') ' 1'
      write(lw,fmt='(A,T111,I6)') ' rm -f zgoubi.OPTICS.out'

      write(lw,fmt='(A,T111,I6)') '''OBJET''',it
      write(lw,FMT='(2(F15.6,9x,A))') 
     >  broDef*1.d3,' ! reference rigidity (kG.cm) ',TGeV,' GeV - kin'
      write(lw,fmt='(A)') '5.1' 
      write(lw,fmt='(A)') '.001 .001 .001 .001 0. .0001  '
C      write(lw,fmt='(A,F16.8)') '0. 0. 0. 0. 0. ',bro
      write(lw,fmt='(A,F16.8)') '0. 0. 0. 0. 0. ', 1.d0
      write(lw,fmt='(1p,4e13.5,A,4e13.5)') 
     >alfx,betx,alfy,bety,' 0. 1. ',dx,dxp,dy,dyp

c      it = it+1
      write(lw,fmt='(A,T111,I6)') '''PARTICUL'''
C      write(lw,fmt='(A)') '9.3827203E+02 1.602176487E-19 1.7928474 0 0 '
      write(lw,fmt='(A)') 
     >       ' 0.51099892 1.60217653e-19 1.15965218076e-3 0. 0.'

c      it = it+1
      write(lw,fmt='(A,T111,I6)') '''SPNTRK'''
      write(lw,fmt='(A)') '4.1 '
      write(lw,fmt='(A)') '0. 0.  1. '

c      it = it+1
      write(lw,fmt='(A)') ' '
      write(lw,fmt='(A,T111,I6)') '''OPTICS'''
      write(lw,fmt='(A,T111,I6)') ' 1   all PRINT'
      write(lw,fmt='(A)') ' '
c      it = it+1
      write(lw,fmt='(A,T111,I6)') '''SRLOSS''' 
      write(lw,fmt='(A)') '0.1         ! srLoss'
      write(lw,fmt='(A)') 'BEND  scale '
      write(lw,fmt='(A)') '1  123456'
c      it = it+1
      write(lw,fmt='(A)') ' '
      write(lw,fmt='(A,T111,I6)') '''SCALING'''
      write(lw,fmt='(A,T111,I6)') ' 1  2'
      write(lw,fmt='(A,T111,I6)') 'MULTIPOL'
      write(lw,fmt='(A)') ' -1'
      write(lw,fmt='(a,F19.9)') ' ',broDef
      write(lw,fmt='(A)') ' 1'
      write(lw,fmt='(A,T111,I6)') 'BEND'
      write(lw,fmt='(A)') ' -1'
      write(lw,fmt='(a,F19.9)') ' ',broDef
      write(lw,fmt='(A)') ' 1'
      write(lw,fmt='(A)') ' '

      return
      end 

      subroutine lmnt
     >(lw,lout,bro,frf,ik,noel,kley,name,xl,ang,ak1,ak2,ak3,ak4,
     >                                        tilt,e1,e2,h1,h2,
     >                                it)
c      subroutine lmnt
c     >(lw,bro,frf,ik,noel,kley,name,xl,ang,ak1,ak2,tilt,e1,e2,h1,h2,
c     >                                                       it)
      implicit double precision (a-h,o-z)
      character kley*(*)

      character name*16, ny*1

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

      integer debstr, finstr

      character txt*80, txtfrm*23
      save kpos, kbm

      logical okV, okVi, Vbnd
      save okV

      pi = 4.d0*atan(1.d0)
      zero=0.d0
      txtfrm = '(A,T12,A,T22,A,T111,i6)'

C----- 1    2    3    4    5    6    7    8    9    10   11   12   13
C----- DRIF RBEN SBEN QUAD SEXT OCTU MULT SROT YROT MARK KICK HMON VMON
C----- 14   15   16   17   18   19   20
C----- HKIC VKIC SOLE RCOL MATR INST MONI
       goto ( 1, 2, 2, 4, 5, 6, 7, 8, 9,
     >       10,11,12,12,11,11,13, 1, 1, 1, 1) ik

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
          write(lw,fmt='(G10.4,G8.2,G16.10,9G7.1)')
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
       endif
C----- R- or SBEN -> MULTIPOL
       if(xl .eq. 0.d0) then
         write(lw,fmt=txtfrm) '''MARKER''',name,kley,it
         goto 99
       endif

!      xlmag = xl
!      ro=xlmag/(2.d0*sin(ang/2.d0))
!      xlarc = ro * ang
      if(ik.eq.2) then
C case RBEN
        xlarc = xl/(2.d0*sin(ang/2.d0)) * ang
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
C          ale = -0.5d0 * ang
          ale = 0.d0
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

      Vbnd = okV .and. tilt.gt.1d-6

      if(Vbnd) then
        it = it+1
        write(lw,fmt=txtfrm) '''TRAROT''',name,'_TRAROT_I',it
        write(lw,fmt='(a,1p,e14.6,a)') '0. 0. 0. ',tilt,' 0. 0.'
c        if(kpos.eq.3 .or. kpos.eq.4) then
c          ale = 0.d0
c          xce=0.d0
c          yce=0.d0
c        endif
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

      if(tilt.gt.1d-6) then
        if(.not. Vbnd) then
          b1 = 0.d0
          b2 = 0.d0
          b3 = 0.d0          
        endif
      endif

      write(lw,fmt='(F9.4,F7.2,3F13.8,7F4.1)')
     >xlmag*cm, x10, b1*10.d0, b2*10.d0, b3*10.d0, x0,x0,x0,x0,x0,x0,x0

      txt=txffd
      if(frf .eq. 0.) txt='0. 0. '//txt
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfd
C      write(lw,fmt='(f14.8,a)') tilt,' 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(a)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      if(xlmag*cm-int(xlmag*cm).gt.0.5) then 
            istepdip=int(xlmag*cm)+1
      else
            istepdip=int(xlmag*cm)
      endif
      if(istepdip.lt.10) then
        write(lw,fmt='(2A)') ' #30|9|30    Dip',name
      elseif(istepdip.ge.10.and.istepdip.lt.100) then
        write(lw,fmt='(A,I2,A,2X,A)') ' #30|',istepdip,'|30    Dip',name
      elseif(istepdip.ge.100.and.istepdip.lt.1000) then
        write(lw,fmt='(A,I3,A,2X,A)') ' #30|',istepdip,'|30    Dip',name
      elseif(istepdip.ge.1000) then
        write(lw,fmt='(2A)') ' #30|9999|30    Dip',name
      endif

      if    ( name(1:2) .eq. 'UD' ) then
        if(b2.gt.0.d0) yce = 1.068472794
        if(b2.lt.0.d0) yce = 1.059533608
      elseif( name(1:2) .eq. 'WD' ) then
        if(b2.gt.0.d0) yce = 1.336991312
        if(b2.lt.0.d0) yce = 1.323022405
      elseif( name(1:2) .eq. 'XD' ) then
        if    (int(xlmag*cm) .eq. 365 ) then
          if(b2.gt.0.d0) yce = -1.480420668
          if(b2.lt.0.d0) yce = -1.455473251
        elseif(int(xlmag*cm) .eq. 294 ) then
          if(b2.gt.0.d0) yce = -0.9577199009
          if(b2.lt.0.d0) yce = -0.9472215129
        endif
      elseif( name(1:2) .eq. 'YD' ) then
        if    (int(xlmag*cm) .eq. 365 ) then
          if(b2.gt.0.d0) yce = 1.480420668
          if(b2.lt.0.d0) yce = 1.455479348
        elseif(int(xlmag*cm) .eq. 294 ) then
          if(b2.gt.0.d0) yce = 0.9577199009
          if(b2.lt.0.d0) yce = 0.9472215129
        endif
      endif

      if(kpos.ne.4) then
        if(tilt.gt.1d-6) then
            write(lw,*)  ' 1 0. 0. 0. '
        else
           write(lw,*)  kpos,xce,yce,ale
        endif
      else
       write(lw,*)  1,0.d0,0.d0,0.d0
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

      if(Vbnd) then
        it = it+1
        write(lw,fmt=txtfrm) '''TRAROT''',name,'_TRAROT_O',it
C        write(lw,fmt=txtfrm) '''TRAROT''',name,kley,it
        write(lw,fmt='(a,1p,e14.6,a)') '0. 0. 0. ',-tilt,' 0. 0.'
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

 4    continue
C----- QUAD
      if(xl .eq. 0.d0) then
        write(lw,fmt=txtfrm) '''MARKER''',name,kley,it
        goto 99
      endif

c         write(*,*) 'QUAD ',ak1,bro,ak1*bro
c              read(*,*)
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,kley,it
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
      if(xl*cm-int(xl*cm).gt.0.5) then 
            istepdip=int(xl*cm)+1
      else
            istepdip=int(xl*cm)
      endif
      istepdip = istepdip
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
        write(lw,fmt=txtfrm) '''MARKER''',name,kley,it
        goto 99
      endif
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .Sext'
      b3 = 10.d0 * ak2 * bro *(x10/100.)**2 / 2.d0
      write(lout,*) ' SEXT         b3 = ',b3,' kG'
      write(lw,fmt='(F10.4,F6.2,2F10.4,F16.10,7F4.1)')
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
        write(lw,fmt=txtfrm) '''MARKER''',name,kley,it
C        it = it - 1
        goto 99
      endif
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .Octu'
      write(lw,fmt='(F10.4,F6.2,2F16.10,4F9.4)')
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
C      if(abs(ak1)+abs(ak2)+abs(ak3)+abs(ak4) .ne. 0.d0) xl = 1.d-5
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
      b3 = ak2 * bro *(x10/100.)**2 / 2.d0  / xl * 10.d0
! octupole field at .1 m (kG)
      b4 = 0.d0
      b5 = 0.d0
      b6 = 0.d0
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .Mult'
      write(lw,fmt='(F10.4,F6.2,2F16.10,8F4.1)')
     >xl*cm,x10,x0,b2,b3,b4,b5,b6,x0,x0,x0,x0
      write(*,*)
     >b2,b3,b4,b5,b6,' b2,b3,b4,b5,b6'
      write(*,*) ak1,' ak1 '
      write(87,*)
     > abs(ak1)+abs(ak2)+abs(ak3)+abs(ak4),' abs(k1)+abs(k2...'
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
        write(lw,fmt=txtfrm) '''TRAROT''',name,'_TRAROT_I',it
C      write(lw,fmt=txtfrm) '''TRAROT''',name,kley,it
      write(lw,fmt='(a,1p,e14.6,a)') ' 0. 0. 0. ',e1,' 0. 0. '
      goto 99

 9    continue
C----- YROT ***
        write(lw,fmt=txtfrm) '''TRAROT''',name,'_TRAROT_I',it
C      write(lw,fmt=txtfrm) '''TRAROT''',name,kley,it
      goto 99

 10   continue
C----- MARK ***
      write(lw,fmt=txtfrm) '''MARKER''',name,kley,it
c      write(*,fmt=txtfrm) '''MARKER''',name,kley,it
c             read(*,*)
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
      write(lw,fmt='(1P,G14.7,G10.4,G14.7)') xl*cm,r0,e1*bro*t2kg
      write(lw,fmt='(2F10.2)') r0, r0
      write(lw,fmt='(A,2X,A)') '1.  Soleno',name
      write(lw,fmt='(A)') '1 0. 0. 0. '
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

      entry lmnkbm(kbmi)
      kbm =kbmi 
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

      IF(LENGTH.EQ.0) RETURN

1     CONTINUE
        DEBSTR=DEBSTR+1
        IF (STRING(DEBSTR:DEBSTR) .EQ. ' ') THEN
          IF(DEBSTR .GE. LENGTH) THEN
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
      CHARACTER(*) STR, STR2
C     ---------------------------------------------------------------
C     .TRUE. if the string STR contains the string STR2 at least once
C     IS = position of first occurence of STR2 in STR 
C     (i.e.,STR(IS:IS+LEN(STR2)-1)=STR2)
C     ---------------------------------------------------------------
      INTEGER DEBSTR,FINSTR
      LNG2 = LEN(STR2(DEBSTR(STR2):FINSTR(STR2)))
      IF(LEN(STR).LT.LNG2 .OR.
     >   (DEBSTR(STR).EQ.0 .AND. FINSTR(STR).EQ.0)
     >     ) GOTO 1
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