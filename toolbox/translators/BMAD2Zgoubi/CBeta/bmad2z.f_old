C23456789012345678901234567890123456789012345678901234567890123456789012
      implicit double precision (a-h,o-z)
      character fnr*80, fnw*80, fnw2*81, txt*132
      character txt590*990,keyword*20,txt20*20,txt40*40
      parameter(lr=9, lw=10, lw2=12, lout=14)
      parameter (mxl=15000)
      dimension iq(mxl)
      character kley*4, name*20, kley1*4
      parameter(mxk=22)
      character kle(mxk)*4, ny*1
      character*35 warn
      character*29 string
      logical TOMANY

      logical ok, strcon
      integer debstr, finstr
      logical okAutoRef
      logical gttext
      character(20) fnwSeg

c      data kle/ 'DRIF', 'RBEN', 'SBEN', 'QUAD', 'SEXT', 'OCTU', 'MULT',
c     >'SROT', 'YROT', 'MARK', 'KICK', 'HMON', 'VMON','HKIC', 'VKIC', 
c     >'SOLE', 'RCOL', 'MATR', 'MONI', 'RFCA', 'MATC', 'PATC' /
      data kle/ 'Drif', 'Rben', 'Sben', 'Quad', 'SEXT', 'OCTU', 'MULT',
     >'SROT', 'YROT', 'Mark', 'KICK', 'HMON', 'VMON','HKIC', 'VKIC', 
     >'SOLE', 'RCOL', 'MATR', 'MONI', 'Rfca', 'Matc', 'Patc' /
      data okAutoRef / .false. /

C      write(*,*) ' name of the survey file to be translated:'
C      read(*,fmt='(A)') fnr
C      write(*,*) ' name of the zgoubi data file:'
C      read(*,fmt='(A)') fnw

      fnw='trad.out'
      fnw2='trad.out2'

      fnr='bmadzg.in'
      if(fnr.ne.'bmadzg.in') then
        write(*,*)
        write(*,*) 'fnr MUST be  ''bmadzg.in'' '
        write(*,*)
      endif
      fnw='trad.out'

      OPEN(UNIT=lout,FILE='bmadzg.out',err=97)
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
            kbm = 2
 7    continue
      if( kbm .ne. 2   ) kbm = 1
      if(kbm.eq.1   ) write(*,*) ' BEND will be used for S/RBEND'
      if(kbm.eq.2) write(*,*) ' MULTIPOLE will be used for S/RBEND'
      call lmnkbm(kbm)

C Get momentum
      read(lr,fmt='(A)',end=85,err=99) txt590
      write(*,*) txt590(1:123)
      read(lr,fmt='(A)',end=85,err=99) txt590
      write(*,*) txt590(1:123)
      read(lr,fmt='(A)',end=85,err=99) txt590
      write(*,*) txt590(1:123)
c      read(lr,fmt='(A)',end=85,err=99) txt590
c      write(*,*) txt590(1:123)
      
 55   continue
      read(lr,fmt='(a)',end=85,err=99) txt590
      read(txt590,*) string
      string = trim(string)
c      write(*,*) ' string = ', string
c      read(*,*)      
      if( string(1:1) .eq. '#') goto 55
      
      read(txt590,*)  ibra, indx, string, string, string, s, t, Etot
        read(txt590(85:),*,end=85,err=95)
     >  sexit, time, Etot, xl, b0, b1, b2, 
     >  rfG, e1, e2, xoff, yoff, zoff, xptch, x1li, x2lim, y1lim, y2lim, 
     >  string, string, betx, alfx, gamx, phix, Dx, Dxp, 
     >                  bety, alfy, gam2, phiy, Dy, Dyp, etaz, 
     >       xco, xpco, yco, ypco
      write(*,*) ' branch  index :', ibra, indx
      backspace(lr)
      
      it = 0
      call objet(lw,it,s,Etot,
     >                           bro,
     >           betx, alfx, Dx, Dxp, 
     >                  bety, alfy, Dy, Dyp, 
     >       xco, xpco, yco, ypco)

      write(*,*) 
      write(*,*) ' Now translating. Busy...'

      ir = 0
      noel = 0
      sxl = 0.d0
      s1 = 0.d0
      s = 0.d0
      bro1 = 1.d0     ! T.m 
      pnlty_indiv = 1.d-12   ! For match of individuel BENDs. This is done on-line. *VERY* imporatnt for zero-ing H and V orbits 
      okAutoref = .false.     ! will insert AUTOREF on particle #1 for use with OBJET/KOBJ=5 or 2. 
      kley1 = ' '
      drft = 0.d0

 86   continue
      
        read(lr,fmt='(a)',end=85,err=99) txt590
        read(txt590,*) string
        string = trim(string)
c        write(*,*) ' string = ', string
c        read(*,*)      
        if( string(1:1) .eq. '#') goto 86

c         write(*,*) txt590(1:123)
c        read(txt590(1:16),*,end=85,err=95)
c     >  int1, int2
c        write(*,*) 'txt590(1:16) : ' ,int1,int2
c        read(txt590(17:88),*,end=85,err=95) name, kley, string
c        read(txt590(89:finstr(txt590)),*,end=85,err=95)
c     >  sexit, time, Etot, xl, b0, b1, b2, 
c     >  rfG, xoff, yoff, zoff, xptch, x1li, x2lim, y1lim, y2lim, 
c     >  string, string, betx, alfx, gamx, phix, Dx, Dxp, 
c     >                  bety, alfy, gam2, phiy, Dy, Dyp, etaz, 
c     >  xco, xpco, yco, ypco, zco, zpco

        read(txt590,*,end=85,err=95)
     >  ibra,indx, name,kley, string
        read(txt590(85:),*,end=85,err=95)
     >  sexit, time, Etot, xl, b0, b1, b2, 
     >  rfG, e1, e2, xoff, yoff, zoff, xptch, x1li, x2lim, y1lim, y2lim, 
     >  string, string, betx, alfx, gamx, phix, Dx, Dxp, 
     >                  bety, alfy, gam2, phiy, Dy, Dyp, etaz, 
     >       xco, xpco, yco, ypco
c , zco, zpco
c sigma_x         sigma_px          sigma_y         sigma_py          sigma_z        
c sigma_pz      norm_emit_x      norm_emit_y      norm_emit_a      norm_emit_b

      write(*,*) ' branch  index :', ibra, indx,sexit
        
        if(kley .eq. 'PIPE') kley = 'DRIF'
        if(kley .eq. 'INST') kley = 'DRIF'
        if(kley .eq. 'MATC') kley = 'MARK'
        if(kley .eq. 'LCAV') kley = 'RFCA'
        if(kley .eq. 'Pipe') kley = 'Drif'
        if(kley .eq. 'Inst') kley = 'Drif'
        if(kley .eq. 'Matc') kley = 'Mark'
        if(kley .eq. 'Lcav') kley = 'Rfca'

        am = 0.510998902E6
        T = Etot -am
        pmom = sqrt(T*(T + 2.d0 * am))
        bro = pmom /2.99792458d8    ! T.m 

        ang = abs(b0 * xl /bro )
        ak1 = b1 / bro

c              write(*,*) ' bmad2z ',name,ak1, b1, bro
c                  read(*,*)

        ak2 = b2 / bro
        ak3 = 0.d0 
        ak4 = 0.d0 
        tilt = 0.d0 
c        e1 = 0.d0 ; e2 = 0.d0
c        h1 = 0.d0 ; h2 = 0.d0 
        sCntr = sexit - xl/2.d0
        if(rfG .ne. 0.d0) volt = rfG * xl 
        fnwSeg = name(debstr(name):debstr(name)+2)        
C RF cavity
        phas = 0.d0 
        freq = 1.3d9
C--

        write(*,*) '-----------------------------------'
        write(*,fmt='(i4,2x,3A,2x,1p,4e12.4)') 
     >  int2,kley, '  ', name, xl, ang, ak1, ak2
              
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
 100    FORMAT(/,10X,' Key ***',A,'*** not translated...')
        goto 86 
C        stop
 
 87     CONTINUE
        it = it + 1
c        call lmnt(
c     >    lw,bro1,frf,ik,noel,kley,name,xl,ang,ak1,ak2,ak3,ak4,
        call lmnt(
     >    lw,bro,frf,ik,noel,kley,name,xl,ang,b0,ak1,ak2,ak3,ak4,
     >                                            tilt,e1,e2,h1,h2,
     >                                      hkic, vkic,
     >                         volt, phas, freq, hgap, fint,pnlty_indiv, 
     >                          okAutoref,  it, fnwSeg, 
     >          xoff, yoff, zoff, xptch)
        ds = s - s1
        s1 = s
        sxl = sxl + xl

c        write(88,fmt='(1p,5(e16.8,1x),2(A16,1x),A35)') 
c     >                         s,sxl,s-sxl,xl,ds,kley,name,warn

      goto 86

 85   continue
      write(lw,fmt='(A,T111,I6)') '! ''TWISS''',it
      write(lw,fmt='(A)')   '! 1  1. 1. '
c      write(lw,fmt='(A,T111,I6)') '''MATRIX''',it
c      write(lw,fmt='(A)')   '1  11 '
c      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''END''',it
C      it = it + 1
      write(lw,fmt='(///,A)') ' '
      write(lw,fmt='(///,A)') '''REBELOTE'''
      write(lw,fmt='(A)')   '600 0.1 99'
      write(*,*) ' Read ',ir,' elements from bmadzg.in file ',fnr
      write(*,*) ' Translated to',it,
     >' elements into the zgoubi data file ', fnw
      write(*,*) ' end of bmadzg.in file '
      goto 999

 95   write(*,*) ' error during read txt590'
      goto 999
 97   write(*,*) ' error open bmadzg.in file'
      goto 999
 99   write(*,*) ' error during read in bmadzg.in file'
      goto 999

 999  continue
      close(lr)
      close(lw)
      close(lw2)
      stop
      end
        
      subroutine objet(lw,it,s,Etot,
     >                            bro,
     >           betx, alfx, Dx, Dxp, 
     >                  bety, alfy, Dy, Dyp, 
     >       xco, xpco, yco, ypco)
      implicit double precision (a-h,o-z)
      parameter (cm=1d2,mrd=1d3)

 1    continue

        am = 0.510998902E6
       T = Etot -am
       pmom = sqrt(T*(T + 2.d0 * am))
        bro = pmom /2.99792458d8    ! T.m 
        gamma = sqrt(pmom*pmom + am*am)/am
      write(*,*) ' Problem rigidity (T.m), momentum (eV/c), G.gamma : ',
     >   bro,pmom,1.79284735d0*gamma
          
      write(lw,*) 'Generated by BMAD -> Zgoubi translator'
      write(lw,fmt='(A,T111,I6)') '''OBJET''',it
      write(lw,FMT='(F19.6,7x,2(2x,A,F17.4))') bro*1.d3,
     >'reference energy (total) = ',Etot, 
     >',  G.gamma = ',1.79284735d0*gamma
      write(lw,fmt='(A)') '5.1' 
      write(lw,fmt='(A)') '.001 .001 .001 .001 0. .0001  '
      write(lw,fmt='(1p,5(e13.5,1x),a)') 
     >xco*cm,xpco*mrd,yco*cm,ypco*mrd,s*1.d2,' 1.'
      write(lw,fmt='(4(f10.5,1x),a,4(f10.4,1x))') 
     >alfx,betx,alfy,bety,'0. 1. ',Dx,Dxp,Dy,Dyp
      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''PARTICUL''',it
      write(lw,FMT='(A)')'0.51099892 1.60217653e-19 0.0011596521811 0 0'

      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''FAISCEAU''',it

c      it = it + 1
c      write(lw,fmt='(A,T111,I6)') '''OPTICS''',it
c      write(lw,FMT='(A)') 
c     > ' 2   Print out transport coeffs to zgoubi_MATRIX_out' 

      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''OPTICS''' ,IT
      write(lw,fmt='(A)') '1 all PRINT'
      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''SCALING''' ,IT
      write(lw,fmt='(A)') '0  6 '
      write(lw,fmt='(A)') 'BEND' 
      write(lw,fmt='(A)') '   -1   ' 
      write(lw,fmt='(1P,E16.8)') BRO
      write(lw,fmt='(A)') '1  '
      write(lw,fmt='(A)') 'MULTIPOL HKIC'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8)') BRO
      write(lw,fmt='(A)') '1      '
      write(lw,fmt='(A)') 'MULTIPOL VKIC'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8)') BRO
      write(lw,fmt='(A)') '1      '
      write(lw,fmt='(A)') 'MULTIPOL QUADQF'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8)') BRO
      write(lw,fmt='(A)') '1      '
      write(lw,fmt='(A)') 'MULTIPOL QUADQD'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8)') BRO
      write(lw,fmt='(A)') '1      '
      write(lw,fmt='(A)') 'MULTIPOL '
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8)') BRO
      write(lw,fmt='(A)') '1      '
      it = it + 1
      write(lw,fmt='(A)') '  ' 
      write(lw,fmt='(A,T111,I6)') '''DRIFT''' , IT
      write(lw,fmt='(A)') ' 0. ' 
      write(lw,fmt='(A)') '  ' 

      return
      end 

      subroutine lmnt(
     >lw,bro,frf,ik,noel,kley,name,xl,ang,b1,ak1,ak2,ak3,ak4,
     >                                        tilt,e1,e2,h1,h2,
     >                                      hkic, vkic,
     >               volt, phas, freq, hgap, fint,pnlty,
     >                       okAutoref,     it,fnwSeg,
     >          xoff, yoff, zoff, xptch)
c      subroutine lmnt
c     >(lw,bro,frf,ik,noel,kley,name,xl,ang,ak1,ak2,tilt,e1,e2,h1,h2,
c     >                                                       it)
      implicit double precision (a-h,o-z)
      character(*) kley, fnwSeg
      logical  okAutoref
      character(*) name
      character ny*1

      character*80 txfd, txfq,fqlhc , fdlhc ,fqrec , fdrec,fd,fq
      character*80 txffd,txffq,ffqlhc,ffdlhc,ffqrec,ffdrec,ffd,ffq
      character*80 ffdmu, fdmu, ffqmua,ffqmus,fqmua
      save txfd,txfq,txffd, txffq
      logical drimul
      save drimul
      character(4) txt4
      character(10) lbl2

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
      save kpos, kbm

      logical okV, okVi, Vbnd
      save okV
      integer debstr, finstr
      logical match
      save match
      character(35) fitxt(50)
      character(6) txt6
      character(2) txt2
      logical ok, idluni
      character(100) cmmnd
      dimension var(50)
      logical gttext
      logical strcon

      data match / .true. /

C     electrons ----------------------------
         b1     = -b1
c        ang  = -ang
c        e1 = -e1
c        e2 = -e2
        ak1=-ak1 ; ak2= -ak2 ; ak3=-ak3 ; ak4=-ak4
C---------------------------------------

      pi = 4.d0*atan(1.d0)
      zero=0.d0
      txtfrm = '(A,T12,A,2X,A,T111,i6)'

C----- 1    2    3    4    5    6    7    8    9    10   11   12   13
C----- DRIF RBEN SBEN QUAD SEXT OCTU MULT SROT YROT MARK KICK HMON VMON
C----- 14   15   16   17   18   12   19   20   22
C----- HKIC VKIC SOLE RCOL MATR MONI RFCA MATC PATC
       goto (
     > 1,   2 ,  2 ,  4,   5,   6,   7,   8,   9,   10,  11,  12,  12,
     > 11,  11,  13,  1,  18,   12,   19, 10, 22) ik

 1    continue
C----- DRIF      
      if(xl .eq. 0.d0) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm)'''MARKER''',name,' '//lbl2,it
        goto 99
      else
        if(drimul) then
          dum=1.d-20
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
          write(lw,fmt=txtfrm) '''MULTIPOL''',name,' '//lbl2,it
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
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
          write(lw,fmt=txtfrm) '''DRIFT''',name,' '//lbl2,it
          write(lw,fmt='(F23.15)') xl*cm
        endif
      endif
      goto 99

 2    continue
       if    (kbm.eq.1) then
C-------- Translate R- or SBEND to BEND
         if(ik.eq.2) then
C---------- RBEN -> BEND
           goto 210
         elseif(ik.eq.3) then
C---------- SBEN -> BEND
           goto 3
         endif
       elseif(kbm.eq.2) then
C----- R- or SBEN -> MULTIPOL
         if(xl .eq. 0.d0) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
           write(lw,fmt=txtfrm) '''MARKER''',name,' '//lbl2,it
           goto 99
         endif

         if(name .eq. 'MXQ1S01') goto 31   !  -> DIPOLE

       elseif(kbm.eq.3) then
C----- R- or SBEN -> DIPOLE
         if(xl .eq. 0.d0) then
           lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
           write(lw,fmt=txtfrm) '''MARKER''',name,' '//lbl2,it
           goto 99
         endif
         goto 31
       endif

          xce = 0.d0
          yce = 0.d0
          ale = 0.d0
          ycea = 0.d0
          yceb = 0.d0
          alea = 0.d0
          aleb = 0.d0

          frf = 0.
C           if(fint*hgap .ne.0.d0) frf =1.          

      if(ik.eq.2) then
C case RBEN
        xlarc = xl/(2.d0*sin(ang/2.d0)) * ang
      elseif(ik.eq.3) then
C case SBEN
        xlarc = xl
      else
        stop ' sbr lmnt, no such option ik '
      endif

      off = 0.d0
      facf = .6700000

      if(ang.ne.0.d0) then
C Ang .ne. 0

        ro = xlarc /ang
        xlmag = 2.d0 * ro * sin (ang/2.d0)
c        b1 = bro / ro                              ! dipole field (T)

        if(abs(e1)+abs(e2) .gt. 1.d-6) then
C     e1 or e2 non-zero

           write(*,*) ' bmad2z e1, e2, ang, b1 ',e1,e2,ang,b1
           
c          if(  abs(e1 - ang) .lt. 1e-6 ) then
          if(  abs(abs(e1) - ang) .lt. 1e-6 ) then
C e1=ang, e2=0
            kpos =1
            xce=0.d0
            ycea= 0.
            yceb= abs( ro*(1.d0-cos(ang))) * cm 
            off = yceb * .6700000 /cm ! chck sign vs. sign(k1)
            off = yceb * facf /cm 
            if(ak1.lt.0.d0) off = -off
C               off=0.
c                 write(*,*) ' yceb =',yceb
c                      stop
            xlmag = ro * sin (ang)
            alea = -abs(e1)
            aleb = 0.

c          elseif(  abs(e2 - ang) .lt. 1e-6 ) then
          elseif(  abs(abs(e2) - ang) .lt. 1e-6 ) then
C e1=0, e2=ang
            kpos =1
            xce=0.d0
            ycea= - abs(ro*(1.d0-cos(ang))) * cm 
            yceb= 0. 
            off = ycea * 0.67000000 /cm ! empirical, for k1>0.  Chck sign vs. sign(k1)
            off = ycea * facf /cm 
C            if(ak1.lt.0.d0) off = -off
            if(ak1.gt.0.d0 .and. b1.lt.0.d0) off = -off
            xlmag = ro * sin (ang)
            alea = 0.
            aleb = -abs(e2)

c            write(*,*) ' abs(abs(e2) - ang) .lt. 1e-6 ',
c     >      abs(abs(e2) - ang) .lt. 1e-6 
c           read(*,*)
            
          else
C e1 = e2
            kpos = 3
            xce=0.d0
            yce=0.d0
            yce=  ro*(1.d0-cos(ang/2)) * cm
            off = yce *.67000000 /cm 
            off = yce *facf /cm 
            off = abs(off)
            if(b1.gt.0. ) off = -off
            yce = 0.
C               off=0.
            ak1 = -ak1
c            if(b1<0.d0) then 
C Case of MXJ1S03 
c              ale = ang/2.d0
c            else
              ale = -ang/2.d0
c            endif
          endif

        else   ! if(abs(e1)+abs(e2) .gt. 1.d-6) then
          if(kpos.eq.4) then
            ale = -0.5d0 * ang * 180/pi
            yce = ro*(1-cos(ale*pi/180)) * cm
          else
            kpos = 3
            xce=0.d0
            yce=0.d0
            yce=  ro*(1.d0-cos(ang/2)) * cm
            off = yce *.67000000 /cm 
            off = yce *facf /cm 
            off = abs(off)
            if(b1.gt.0. ) off = -off
            yce = 0.
C               off=0.
            ak1 = -ak1
c            if(b1<0.d0) then 
C Case of MXJ1S03 
c              ale = ang/2.d0
c            else
              ale = -ang/2.d0
c            endif
          endif
        endif

      else
C Ang = 0
        xlmag=xl
        b1=0.d0
          ale=0.d0
          xce=0.d0
          yce=0.d0
      endif

      Vbnd = okV .and. abs(tilt).gt.1d-6

c             write(*,*) ' 2 ',name, ang, b1, kpos, Vbnd

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

          b1mod = b1 + ak1 * off * Bro    ! Tesla

c             write(*,*) ' b1mod ',b1,b1mod


      endif
! quadrupole field at .1 m (T)
      b2 = ak1*bro * x10/100.                     
! sextupole field at .1 m (T)
      b3 = ak2 * bro *(x10/100.)**2 / 2.d0
C         write(*,*) ak2, bro, b3, x10

      if(abs(tilt).gt.1d-6) then
        if(.not. Vbnd) then
c          b1mod = 0.d0
          b2 = 0.d0
          b3 = 0.d0          
        endif
      endif
      
      absb1 = abs(b1)

      xes = 2.5 * gap

      YCE3 = 0.d0
      bdip =  absb1*10.d0
      dipk = b2*10.d0

      if(hgap .ne. 0.d0) then
        gap = 2.*hgap
        if (gap .lt. 0.08) gap = 0.08
        if(3.*gap .gt. xlmag/2.) then
          gap = xlmag/2.1 /3.
        endif 
      endif
      FFdp = gap*cm

           match = .false.

      if(match) then

        call system('mkdir -p ./bendMatch/')
        ok = idluni(lunW)
        open(unit=lunW,file=
     >  './bendMatch/'//name(debstr(name):finstr(name))//'_FIT.dat')
        write(lunw,*) 'Generated using elegant -> Zgoubi translator'
        write(lunw,fmt='(A,T111,I6)') '''OBJET'''
        write(lunw,FMT='(a)') ' 1000. '
        write(lunw,fmt='(A)') '5' 
        write(lunw,fmt='(A)') '.01 .01 .01 .01 0. .001  '
        write(lunw,fmt='(A)') '0. 0. 0. 0. 0. 1. '
        nuel = 1
        nv = 0

        if(Vbnd) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
          write(lunw,fmt=txtfrm) '''TRAROT''',name,' '//lbl2,it
          write(lunw,fmt='(a,1p,e16.8,a)') '0. 0. 0. ',tilt,' 0. 0.'
          nuel = nuel + 1
        endif

        if(b1<0.d0) then
          write(lunw,fmt=txtfrm) '''YMY''','YMY_IN ','NEG_B',it
          nuel = nuel + 1
        endif
        
        if    (kpos.eq.1) then
          write(lunw,fmt=txtfrm) '''CHANGREF''',' ',' ',it
          write(lunw,fmt='(1p,2(1x,A,e14.6))') 
     >           'YS ',ycea,'ZR ',alea/pi*180.d0
          nuel = nuel + 1
          if(ycea.ne.0.) then
            nv = nv + 1
            write(txt6,fmt='(i0)') nuel
            write(txt2,fmt='(i0)') nv
            fitxt(nv) = txt6//' 1  0  .5      ! var#'//txt2//' ycea'
          endif
        endif

        if(b1<0.d0) then
          write(lunw,*) '''MULTIPOL'''//' '//name//' '//' NEG_B'
        else
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
            write(lunw,fmt=txtfrm) '''MULTIPOL''',name,' '//lbl2,it
        endif

c        write(*,fmt=txtfrm) ' '
c        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
c        write(*,fmt=txtfrm) '''MULTIPOL''',name,' '//lbl2,it
        write(lunw,fmt='(I1,A)') i0,'  .Dip'

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

          b1mod = b1 + ak1 * off * Bro    ! Tesla
c             write(*,*) '  2 b1mod ',b1,b1mod

        else
          b1mod = b1
        endif
! quadrupole field at .1 m (T)
        b2 = ak1*bro * x10/100.                     
! sextupole field at .1 m (T)
        b3 = ak2 * bro *(x10/100.)**2 / 2.d0

        if(abs(tilt).gt.1d-6) then
         if(.not. Vbnd) then
          b1mod = 0.d0
          b2 = 0.d0
          b3 = 0.d0          
         endif
        endif
      
        absb1 = abs(b1mod)

        write(lunw,fmt='(F11.6,F7.2,3F13.8,7F4.1)')
     >  xlmag*cm, x10, absb1*10.d0, b2*10.d0, 
     >  b3*10.d0, x0,x0,x0,x0,x0,x0,x0

        nuel = nuel + 1
        nv = nv + 1
        write(txt6,fmt='(i0)') nuel
        write(txt2,fmt='(i0)') nv
        if(b1.ne. 0.d0) then
          fitxt(nv) = txt6//' 4  0  .5      ! var#'//txt2//' bdip'
c             write(*,*) ' b1, fitxt : ',b1mod,fitxt(nv),nv
        else
          fitxt(nv) = txt6//' 4 0 [-1.,1.]  ! var#'//txt2//' bdip'
c               write(*,*) ' b1, fitxt : ',b1mod,fitxt(nv),nv
        endif

c             write(*,*) name, b1,fitxt(nv)

        nv = nv + 1
        write(txt6,fmt='(i0)') nuel
        write(txt2,fmt='(i0)') nv
        if(b2.ne. 0.d0) then
          fitxt(nv) = txt6//' 5 0 [-3.,3.]  ! var#'//txt2//' kdip'
c               write(*,*) ' b2, fitxt : ',b2,fitxt(nv),nv
        else
          fitxt(nv) = txt6//' 5 0 [-3.,3.]  ! var#'//txt2//' kdip'
c               write(*,*) ' b2, fitxt : ',b2,fitxt(nv),nv
        endif
CCCCCCCCCCCCCCCCCCCCCCc      endif

        if(hgap .ne. 0.d0) then
          gap = 2.*hgap
          if (gap .lt. 0.08) gap = 0.08
          if(3.*gap .gt. xlmag/2.) then
            gap = xlmag/2.1 /3.
          endif 
            nv = nv+1
            write(txt6,fmt='(i0)') nuel
            write(txt2,fmt='(i0)') nv
            fitxt(nv) = txt6//' 15 '//txt6(debstr(txt6):finstr(txt6))
     >          //'.033 .2 ! var#'//txt2//' FFdp'
        endif
        xes = 2.5 * gap
        write(txt,fmt='(f6.2,1x,f6.2)') xes*cm,gap*cm
        txt = txt(debstr(txt):finstr(txt))//
     >  ' 1.00 0.00 0.00 0.00 0.00 0. 0. 0. 0.'
         write(lunw,fmt='(A)') txt(debstr(txt):finstr(txt))
          write(lunw,fmt='(A)') txfd
           write(lunw,fmt='(A)') txt(debstr(txt):finstr(txt))
     >  //' '//name//' FFout'
          write(lunw,fmt='(A)') txfd
         write(lunw,fmt='(a)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
        if(xlmag*cm-int(xlmag*cm).gt.0.5) then 
            istepdip=int(xlmag*cm)+1
        else
            istepdip=int(xlmag*cm)
        endif

        write(lunw,fmt='(A,I0,A)') ' #30|',istepdip,'|30    Dip '//name

c             write(*,*) ' 2 ',name, ang, b1, kpos

        if(kpos .eq. 1) then
                write(lunw,*)  ' 1 0. 0. 0. '
        elseif(kpos.ne.4) then
          if(abs(tilt).gt.1d-6) then
c             write(lunw,*)  ' 1 0. 0. 0. '
             write(lunw,*)  kpos,xce,yce,ale
          else
               write(lunw,*)  kpos,xce,yce,ale
          endif
          if(kpos.eq.3) then
              nv = nv+1
              write(txt6,fmt='(i0)') nuel
              write(txt2,fmt='(i0)') nv
              fitxt(nv) = txt6//' 65 0 [-3.,3]  ! var#'//txt2//' YCE3'
          endif
        else
C             write(lunw,*)  1,0.d0,0.d0,0.d0
             write(lunw,*)  kpos,xce,yce,ale
        endif

        if    (kpos.eq.1) then
            write(lunw,fmt=txtfrm) '''CHANGREF''',' ',' ',it
              write(lunw,fmt='(1p,2(A,e14.6))') 
     >              ' YS ',yceb,' ZR ',aleb/pi*180.d0
          nuel = nuel + 1
          if(yceb.ne.0.) then
            nv = nv+1
            write(txt6,fmt='(i0)') nuel
            write(txt2,fmt='(i0)') nv
            fitxt(nv) = txt6//' 1  0  .5      ! var#'//txt2//' yceb'
          endif
        endif

        if(b1<0.d0) then
             write(lunw,fmt=txtfrm) '''YMY''','YMY_OUT',' ',it
         endif
      
        if(Vbnd) then
          lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
             write(lunw,fmt=txtfrm) '''TRAROT''',name,' '//lbl2,it
              write(lunw,fmt='(a,1p,e16.8,a)') 
     >           '0. 0. 0. ',-tilt,' 0. 0.'
        endif

        write(lunw,FMT='(A)') '''MATRIX'''
        write(lunw,FMT='(A)') ' 1 0 '
        write(lunw,FMT='(A)') '''END'''
        close(lunW)

        cmmnd = '(cd ./bendMatch ; \cp '//
     >  name(debstr(name):finstr(name))//'_FIT.dat zgoubi.dat ; '
     >  //'  ~/zgoubi/current/zgoubi/zgoubi )'
        call system(cmmnd)
        cmmnd = '(cd ./bendMatch ; \cp zgoubi.FITVALS.out '//
     >  name(debstr(name):finstr(name))//'_zgoubi.FITVALS.out )'
        call system(cmmnd)
        cmmnd = '(cd ./bendMatch ; mv -f zgoubi.res '//
     >  name(debstr(name):finstr(name))//'_FIT.res )'
        call system(cmmnd)

C Now get new parameter values as of FIT
        ok = idluni(
     >              lrfit)
        open(unit=lrfit,file='./bendMatch'//'/'//
     >  name(debstr(name):finstr(name))//'_zgoubi.FITVALS.out')
        ok = gttext(6,lrfit,'STATUS OF V',
     >                               txt)
        read(lrfit,*) txt
c        write(*,*) ' gttext  ',txt
c              read(*,*)        
     
        read(lrfit,fmt='(a)',end=444,err=444) txt
c        write(*,*) ' variables  ',txt
c        write(*,*) ' variables  ',txt
c        write(*,*) ' variables  ',txt
        do while (txt(debstr(txt):debstr(txt)+5) .ne. 'STATUS')
c          write(*,*) ' variables  ',txt(debstr(txt):debstr(txt)+5)
          read(txt,*) nlm,nvar,mpar,dum,dum,var(nvar)
c          write(*,*) ' '
c          write(*,*) ' variable value : ', nlm,nvar,mpar,var(nvar)
          read(lrfit,fmt='(a)',end=444,err=444) txt
c          write(*,*) ' variables loop ',txt
c          read(*,*)
        enddo
 444    continue
        close(lrfit)

        do i = 1, nv
           txt4=fitxt(i)(finstr(fitxt(i))-3:finstr(fitxt(i)))
           if(txt4 .eq. 'ycea') ycea=var(i)   
           if(txt4 .eq. 'bdip') bdip=var(i)
           if(txt4 .eq. 'kdip') dipk=var(i)
           if(txt4 .eq. 'FFdp') FFdp=var(i)
           if(txt4 .eq. 'YCE3') YCE3=var(i)   
           if(txt4 .eq. 'yceb') yceb=var(i)   
c            write(*,*) fitxt(i)(finstr(fitxt(i))-3:finstr(fitxt(i)))
c                 write(*,*) ' i, var ', i,var(i) ,name,txt4
c                 write(*,*) ycea,bdip,dipk,FFdp,YCE3,yceb,txt4
        enddo

C                  read(*,*)
C-------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCc      
      endif ! match

      if(Vbnd) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm) '''TRAROT''',name,' '//lbl2,it
        write(lw,fmt='(a,1p,e16.8,a)') '0. 0. 0. ',tilt,' 0. 0.'
          it = it+1
      endif

      if(b1<0.d0) then
        write(lw,fmt=txtfrm) '''YMY''','YMY_IN ','NEG_B',it
          it = it+1
      endif

      if    (kpos.eq.1) then
        if(abs(ycea) + abs(alea) .gt. 1.d-6) then
          write(lw,fmt=txtfrm) '''CHANGREF''',' ',' ',it
          write(lw,fmt='(1p,2(1x,A,e14.6))') 
     >           'YS ',ycea,'ZR ',alea/pi*180.d0
          it = it+1
        endif
      endif

      if(b1<0.d0) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm) '''MULTIPOL''',name,' '//lbl2//' NEG_B',it
      else
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm) '''MULTIPOL''',name,' '//lbl2,it
      endif

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

      endif
! quadrupole field at .1 m (T)
      b2 = ak1*bro * x10/100.                     
! sextupole field at .1 m (T)
      b3 = ak2 * bro *(x10/100.)**2 / 2.d0

      if(abs(tilt).gt.1d-6) then
        if(.not. Vbnd) then
          b2 = 0.d0
          b3 = 0.d0          
        endif
      endif
      
      write(lw,fmt='(F11.6,F7.2,3F13.8,7F4.1)')
     >xlmag*cm, x10, bdip, dipk, 
     >b3*10.d0, x0,x0,x0,x0,x0,x0,x0
      xes = 2.5 * gap
C      write(txt,fmt='(f6.2,1x,f6.2)') xes*cm,gap*cm
      write(txt,fmt='(f6.2,1x,f8.4)') xes*cm,FFdp
      txt = txt(debstr(txt):finstr(txt))//
     >' 1.00 0.00 0.00 0.00 0.00 0. 0. 0. 0.'
      write(lw,fmt='(A)') txt(debstr(txt):finstr(txt))
     >//' '//name//' FFin'
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(A)') txt(debstr(txt):finstr(txt))
     >//' '//name//' FFout'
      write(lw,fmt='(A)') txfd
      write(lw,fmt='(a)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      if(xlmag*cm-int(xlmag*cm).gt.0.5) then 
            istepdip=int(xlmag*cm)+1
      else
            istepdip=int(xlmag*cm)
      endif

      write(lw,fmt='(A,I0,A)') ' #30|',istepdip,'|30    Dip '//name

      if(kpos .eq. 1) then
                write(lw,*)  ' 1 0. 0. 0. '
      elseif(kpos.ne.4) then
        if(abs(tilt).gt.1d-6) then
             write(lw,*)  kpos,xce,YCE3,ale
c             write(lw,*)  ' 1 0. 0. 0. '
        else
               write(lw,*)  kpos,xce,YCE3,ale
        endif
      else
c             write(lw,*)  1,0.d0,0.d0,0.d0
               write(lw,*)  kpos,xce,yce,ale
      endif



c        write(*,*) ' kpos = ',kpos
c        read(*,*)

      if    (kpos.eq.1) then
        if(abs(yceb) + abs(aleb) .gt. 1.d-6) then
          it = it+1
          write(lw,fmt=txtfrm) '''CHANGREF''',' ',' ',it
          write(lw,fmt='(1p,2(A,e14.6))') 
     >              ' YS ',yceb,' ZR ',aleb/pi*180.d0
        endif
      endif

      if(b1<0.d0) then
        it = it + 1
        write(lw,fmt=txtfrm) '''YMY''','YMY_OUT',' ',it
       endif
      
      if(Vbnd) then
        it = it+1
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm) '''TRAROT''',name,' '//lbl2,it
        write(lw,fmt='(a,1p,e16.8,a)') '0. 0. 0. ',-tilt,' 0. 0.'
      endif

      goto 99
 
 210   continue
C----- RBEN -> BEND
      if(xl * ang .eq. 0.d0) then
        if(xl .eq. 0.d0) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
          write(lw,fmt=txtfrm) '''MARKER''',name,' '//lbl2,it
        else
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
          write(lw,fmt=txtfrm) '''DRIFT''',name,' '//lbl2,it
          write(lw,fmt='(1P,E14.6)') xl*cm
        endif
        goto 99
      endif

      Vbnd = okV .and. abs(tilt).gt.1d-6
      if(Vbnd) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm) '''TRAROT''',name,' '//lbl2,it
        it = it+1
        write(lw,fmt='(a,1p,e16.8,a)') '0. 0. 0. ',tilt,' 0. 0.'
      endif

        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
      write(lw,fmt=txtfrm) '''BEND''',name,' '//lbl2,it
      write(lw,fmt='(I1,A)') i0,'  .Bend'
      ro=xl/2./sin(ang/2.)
      b=ang * bro / xl *10.  !bro/ro*t2kg
C      xxl=ro*ang*cm
      ro=ro*cm

C err. FM June 2009      write(lw,fmt='(2F14.7,F15.8)') xl,tilt,b 
      write(lw,fmt='(2F16.9,F15.8)') xl*cm,tilt,b 
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
        write(lw,fmt='(A)') '3 0. 0. 0.'
      else
        write(lw,fmt='(1P,I2,3(1x,e18.10))') kpos,xce,yce,ale
      endif

      if(Vbnd) then
        it = it+1
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm) '''TRAROT''',name,' '//lbl2,it
        write(lw,fmt='(a,1p,e16.8,a)') '0. 0. 0. ',-tilt,' 0. 0.'
      endif

      goto 99

 3    continue
C----- SBEN -> BEND
      if(xl * ang .eq. 0.d0) then
        if(xl .eq. 0.d0) then
          lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
          write(lw,fmt=txtfrm) '''MARKER''',name,' '//lbl2,it
        else
          lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
          write(lw,fmt=txtfrm) '''DRIFT''',name,' '//lbl2,it
          write(lw,fmt='(1P,E14.6)') xl*cm
        endif
        goto 99
      endif

      Vbnd = okV .and. abs(tilt).gt.1d-6
      if(Vbnd) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm) '''TRAROT''',name,' '//lbl2,it
        it = it+1
        write(lw,fmt='(a,1p,e16.8,a)') '0. 0. 0. ',tilt,' 0. 0.'
      endif

cC test
      ro=xl/ang
      b=bro/ro*t2kg

c        ro = xlarc /ang
        xlmag = 2.d0 * ro * sin (ang/2.d0)
c        b1 = bro / ro                              ! dipole field (T)

c            write(*,*) trim(name), e1,e2,ang
c                 read(*,*)
        
        if(abs(e1)+abs(e2) .gt. 1.d-6) then
C     e1 or e2 non-zero

          if(  abs(e1 - ang) .lt. 1e-6 ) then
C e1=ang, e2=0
            kpos =1
            xce=0.d0
            ycea= 0.
            yceb= abs( ro*(1.d0-cos(ang))) * cm 
            off = yceb * .6700000 /cm ! chck sign vs. sign(k1)
            off = yceb * facf /cm 
            if(ak1.lt.0.d0) off = -off
C               off=0.
c                 write(*,*) ' yceb =',yceb
c                      stop
            xlmag = ro * sin (ang)
            alea = -abs(e1)
            aleb = 0.
            te = ang
            ts = 0.
            
          elseif(  abs(e2 - ang) .lt. 1e-6 ) then
C e1=0, e2=ang
            kpos =1
            xce=0.d0
            ycea= - abs(ro*(1.d0-cos(ang))) * cm 
            yceb= 0. 
            off = ycea * 0.67000000 /cm ! empirical, for k1>0.  Chck sign vs. sign(k1)
            off = ycea * facf /cm 
C            if(ak1.lt.0.d0) off = -off
            if(ak1.gt.0.d0 .and. b.lt.0.d0) off = -off
            xlmag = ro * sin (ang)
            alea = 0.
            aleb = -abs(e2)
            te=0
            ts=ang
            
          else
C e1 = e2
            kpos = 3
            xce=0.d0
            yce=0.d0
            yce=  ro*(1.d0-cos(ang/2)) * cm
            off = yce *.67000000 /cm 
            off = yce *facf /cm 
            off = abs(off)
            if(b.gt.0. ) off = -off
            yce = 0.
C               off=0.
            ak1 = -ak1
            if(b<0.d0) then 
C Case of MXJ1S03 
              ale = ang/2.d0
            else
              ale = -ang/2.d0
            endif
          endif

        else
          if(kpos.eq.4) then
            ale = -0.5d0 * ang * 180/pi
            yce = ro*(1-cos(ale*pi/180)) * cm
          else
            kpos = 3
            xce=0.d0
            yce=0.d0
            yce=  ro*(1.d0-cos(ang/2)) * cm
            off = yce *.67000000 /cm 
            off = yce *facf /cm 
            off = abs(off)
            if(b.gt.0. ) off = -off
            yce = 0.d0
C               off=0.
            ak1 = -ak1
            if(b<0.d0) then 
C Case of MXJ1S03 
              ale = ang/2.d0
            else
              ale = -ang/2.d0
            endif
            te=ang/2.
            ts=ang/2.
          endif
        endif

      if(b<0.d0) then
        write(lw,fmt=txtfrm) '''YMY''','YMY_IN ','NEG_B',it
        it = it+1
c        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
c        write(lw,fmt=txtfrm) '''BEND''',fnwSeg,'NEG_B',it
      else
c        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
c        write(lw,fmt=txtfrm) '''BEND''',name,' '//lbl2,it
      endif
      lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
      write(lw,fmt=txtfrm) '''BEND''',name,' '//lbl2,it
      write(lw,fmt='(I1,A)') i0,'  .Bend'
c      xxl = 2.d0*ro*sin(xl/2.d0/ro)
      write(lw,fmt='(F14.7,a,F15.8)') xlmag*cm,' 0. ',abs(b) 
      if(frf .eq. 0.) then
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,te
      else
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,te
      endif
      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
      if(frf .eq. 0.) then
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,ts
      else
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,ts
      endif
      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
C      istepdip=int(xxl*cm/3.0d0)
        write(lw,fmt='(2A)') '0.2    Bend',name 
      if(kpos.eq.3) then
C        write(lw,fmt='(A)') '3 0. 0. 0.'
        write(lw,fmt='(A,1p,e16.8)') '3 0. 0. ',ale
      else
        write(lw,fmt='(1P,I1,2X,3G18.10)') kpos,xce,yce,ale
      endif

      if(b<0.d0) then
        it = it+1
        write(lw,fmt=txtfrm) '''YMY''','YMY_OUT',' ',it
      endif
      if(Vbnd) then
        it = it+1
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm) '''TRAROT''',name,' '//lbl2,it
        write(lw,fmt='(a,1p,e16.8,a)') '0. 0. 0. ',-tilt,' 0. 0.'
      endif

      goto 99

 31   continue
C----- SBEN -> DIPOLE
      if(xl * ang .eq. 0.d0) then
        if(xl .eq. 0.d0) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
          write(lw,fmt=txtfrm) '''MARKER''',name,' '//lbl2,it
        else
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
          write(lw,fmt=txtfrm) '''DRIFT''',name,' '//lbl2,it
          write(lw,fmt='(1P,E14.6)') xl*cm
        endif
        goto 99
      endif

      b1=bro/(xl/ang)*t2kg
      ro=abs(xl/ang)
      ffang = 4.d0     ! deg
      at = abs(ang*180/(4.d0*atan(1.d0))) + 2.*ffang
      acn = ffang
      omgp = acn - ffang
      omgm = - ( at -acn -ffang)
      rm = ro*1.d2   ! cm
      b2 = ro**2 * ak1
      xs = -rm*tan(ACN/180.d0*(4.d0*atan(1.d0)))


      if(match) then

        call system('mkdir -p ./bendMatch/')
        ok = idluni(lunW)
        open(unit=lunW,file=
     >  './bendMatch/'//name(debstr(name):finstr(name))//'_FIT.dat')
        write(lunw,*) 'Generated using elegant -> Zgoubi translator'
        write(lunw,fmt='(A,T111,I6)') '''OBJET'''
        write(lunw,FMT='(a)') ' 1000. '
        write(lunw,fmt='(A)') '5' 
        write(lunw,fmt='(A)') '.01 .01 .01 .01 0. .001  '
        write(lunw,fmt='(A)') '0. 0. 0. 0. 0. 1. '

        nuel = 1
        nv = 0

        Vbnd = okV .and. abs(tilt).gt.1d-6
        if(Vbnd) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
          write(lunw,fmt=txtfrm) '''TRAROT''',name,' '//lbl2,it
          write(lunw,fmt='(a,1p,e16.8,a)') '0. 0. 0. ',tilt,' 0. 0.'
          nuel = nuel + 1
        endif

        if(b1<0.d0) then
          hnorm = -b1
          write(lunw,fmt=txtfrm) '''YMY''','YMY_IN ','NEG_B',it
          nuel = nuel + 1
        else
          hnorm = b1
        endif

c        it = it + 1
        nuel = nuel + 1
        write(lunw,fmt=txtfrm) '''CHANGREF''',name,'      ',it
        write(lunw,fmt='(a,f14.7,t65,a)') ' XS ',xs,
     >' ! RM*tan(ACN-omga+ ==ACN)'   
c        nuel = nuel + 1
c        nv = nv + 1
c        write(txt6,fmt='(i0)') nuel
c        write(txt2,fmt='(i0)') nv
c        if(xs.ne.0.) then
c          fitxt(nv) = txt6//' 1  0  .5      ! var#'//txt2//' XS'
c        else
c          fitxt(nv) = txt6//' 1  0  [-20,20]   ! var#'//txt2//' XS'
c        endif

        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lunw,fmt=txtfrm) '''DIPOLE''',name,' '//lbl2,it
        nuel = nuel + 1
        write(lunw,fmt='(I1,A)') i0,'  .Dipole'
        write(lunw,fmt='(F14.7,2x,F15.8,T65,A)')  AT, RM,' ! AT, RM'
        write(lunw,fmt='(F14.7,2x,2F15.8,a,T65,a)')acn,hnorm,b2,' 0. 0.'
     >  ,'        ! ACNT,  HNORM, indices'
        nv = nv + 1
        write(txt6,fmt='(i0)') nuel
        write(txt2,fmt='(i0)') nv
        fitxt(nv) = txt6//' 5  0  .1      ! var#'//txt2//' hnrm'        
        nv = nv + 1
        write(txt6,fmt='(i0)') nuel
        write(txt2,fmt='(i0)') nv
        if(b2.ne.0.d0) then
          fitxt(nv) = txt6//' 6  0  .9      ! var#'//txt2//' Indx'        
        else
          fitxt(nv) = txt6//' 6  0 [-.1,.1] ! var#'//txt2//' Indx'        
        endif
        write(lunw,fmt='(a)') ' 8. -1.         face 1 ' 
        write(lunw,fmt='(a)') ' 4 .1455 2.2670 -.6395 1.1558 0. 0. 0.'
        e1deg = e1/(4.*atan(1.d0))*180.d0
        write(lunw,fmt='(2(F14.7,2x),a)')omgp,-e1deg,
     >  ' 1.E6 -1.E6 1.E6 1.E6'
        write(lunw,fmt='(a)') ' 8. -1.         face 2 ' 
        write(lunw,fmt='(a)') ' 4 .1455 2.2670 -.6395 1.1558 0. 0. 0.'
        e2deg = e2/(4.*atan(1.d0))*180.d0
        write(lunw,fmt='(2(F14.7,2x),a)')omgm,-e2deg,
     >  ' 1.E6 -1.E6 1.E6 1.E6'
        write(lunw,fmt='(a)') ' 0.  0.         face 3 ' 
        write(lunw,fmt='(a)') ' 4 .1455 2.2670 -.6395 1.1558 0. 0. 0.'
        write(lunw,fmt='(a)') ' 0. 0.   1.E6 -1.E6 1.E6 1.E6  0.'
        write(lunw,fmt='(2A)') ' 2   64           '
        write(lunw,fmt='(2A)') ' .2               '
        rers = rm/cos(ffrad)
        ffrad = ffang/180.*(4.*atan(1.d0))
        write(lunw,fmt='(a,4(f13.7,1x),t70,A)') ' 2  '
     >  ,rers, -ffrad, rers, +ffrad, '  ! 63-67'
        nv = nv + 1
        write(txt6,fmt='(i0)') nuel
        write(txt2,fmt='(i0)') nv
        fitxt(nv) = txt6//' 64 4.066  .2  ! var#'//txt2//' ReRs'        

c        it = it + 1
        write(lunw,fmt=txtfrm) '''CHANGREF''',name,'      ',it
        write(lunw,fmt='(a,f14.7,t65,a)') ' XS ',xs,
     >  ' ! RM*tan(ACN-omga+ ==ACN)'   
        nuel = nuel + 1

        if(b<0.d0) then
          write(lunw,fmt=txtfrm) '''YMY''','YMY_OUT',' ',it
        endif

        if(Vbnd) then
c          it = it+1
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
          write(lunw,fmt=txtfrm) '''TRAROT''',name,' '//lbl2,it
          write(lunw,fmt='(a,1p,e16.8,a)') '0. 0. 0. ',-tilt,' 0. 0.'
        endif

        write(lunw,FMT='(A)') '''MATRIX'''
        write(lunw,FMT='(A)') ' 1 0 '
        write(lunw,fmt='(A)') '''FIT2'' '
        write(lunw,FMT='(i0,a)') nv , '    noSYSout save'
        do i = 1, nv
          write(lunW,*) fitxt(i)
        enddo

        if(len(name(debstr(name):finstr(name))).gt.1)
     >   call getRijInd(name(debstr(name):finstr(name)),
     >                   r11,r12,r33,r34)
        nc = 8
        write(lunw,FMT='(i0,1x,1p,e10.2)') nc,pnlty
        write(lunw,FMT='(a,1p,e13.5,a)') '1 1 1 #End ',R11,' 1. 0'
        write(lunw,FMT='(a,1p,e13.5,a)') '1 1 2 #End ',R12,' 1. 0'
        write(lunw,FMT='(a,1p,e13.5,a)') '1 3 3 #End ',R33,' 1. 0'
        write(lunw,FMT='(a,1p,e13.5,a)') '1 3 4 #End ',R34,' 1. 0'
        write(lunw,FMT='(a)') '3 1 2 #End 0. 1. 0'
        write(lunw,FMT='(a)') '3 1 3 #End 0. 1. 0'
        write(lunw,FMT='(a)') '3 1 4 #End 0. 1. 0'
        write(lunw,FMT='(a)') '3 1 5 #End 0. 1. 0'
        write(lunw,FMT='(A)') '''MATRIX'''
        write(lunw,FMT='(A)') ' 1 0 '
        write(lunw,FMT='(A)') '''END'''
        close(lunW)

        cmmnd = '(cd ./bendMatch ; \cp '//
     >  name(debstr(name):finstr(name))//'_FIT.dat zgoubi.dat ; '
     >//'  ~/zgoubi/current/zgoubi/zgoubi )'
        call system(cmmnd)
        cmmnd = '(cd ./bendMatch ; \cp zgoubi.FITVALS.out '//
     >  name(debstr(name):finstr(name))//'_zgoubi.FITVALS.out )'
        call system(cmmnd)
        cmmnd = '(cd ./bendMatch ; mv -f zgoubi.res '//
     >  name(debstr(name):finstr(name))//'_FIT.res )'
        call system(cmmnd)

C Now get new parameter values as of FIT
        ok = idluni(
     >              lrfit)
        open(unit=lrfit,file='./bendMatch'//'/'//
     >  name(debstr(name):finstr(name))//'_zgoubi.FITVALS.out')
        ok = gttext(6,lrfit,'STATUS OF V',
     >                               txt)
        read(lrfit,*) txt
c        write(*,*) ' gttext  ',txt
c              read(*,*)        
     
        read(lrfit,fmt='(a)',end=445,err=445) txt
c        write(*,*) ' variables  ',txt
c        write(*,*) ' variables  ',txt
c        write(*,*) ' variables  ',txt
        do while (txt(debstr(txt):debstr(txt)+5) .ne. 'STATUS')
c          write(*,*) ' variables  ',txt(debstr(txt):debstr(txt)+5)
          read(txt,*) nlm,nvar,mpar,dum,dum,var(nvar)
c          write(*,*) ' '
c          write(*,*) ' variable value : ', nlm,nvar,mpar,var(nvar)
          read(lrfit,fmt='(a)',end=445,err=445) txt
c          write(*,*) ' variables loop ',txt
c          read(*,*)
        enddo
 445    continue
        close(lrfit)

        do i = 1, nv
           txt4=fitxt(i)(finstr(fitxt(i))-3:finstr(fitxt(i)))
           if(txt4 .eq. 'hnrm') hnorm =var(i)   
           if(txt4 .eq. 'Indx') b2 =var(i)
           if(txt4 .eq. 'ReRs') ReRs=var(i)

c           write(*,*) ' e2z i, txt4, var(i) : ',i, txt4,hnorm,b2,ReRs
c                   read(*,*)

        enddo

C-------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCc      
      endif ! match

      Vbnd = okV .and. abs(tilt).gt.1d-6
      if(Vbnd) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm) '''TRAROT''',name,' '//lbl2,it
        it = it+1
        write(lw,fmt='(a,1p,e16.8,a)') '0. 0. 0. ',tilt,' 0. 0.'
      endif

      if(b1<0.d0) then
c        hnorm = -b1
        write(lw,fmt=txtfrm) '''YMY''','YMY_IN ','NEG_B',it
          it = it+1
c      else
c        hnorm = b1
      endif

c         write(*,*) ' e2z i, txt4, var(i) //// ',i, txt4,hnorm,b2,ReRs
c                   read(*,*)

      write(lw,fmt=txtfrm) '''CHANGREF''',name,'      ',it
      it = it + 1
      write(lw,fmt='(a,f14.7,t65,a)') ' XS ',xs,
     >' ! RM*tan(ACN-omga+ ==ACN)'   

c 'TRAROT'   MXQ1S01               SBEN             
c 0. 0. 0.  -1.570796E+00 0. 0.
c 'CHANGREF'
c XS -22.0023076169     ! RM*tan(ACN-omga+ ==ACN)
c 'DIPOLE' MXQ1S01               SBEN                
c  0                                                                             
c 26.6341335016  314.647658107                             AT, RM                               
c 4. 3.177570990 2.9806237063E-02 0. 0.                          ACNT,  HNORM, indices           
c 8. -1.                                       face 1                    ! 10                      
c 4  .1455   2.2670  -.6395  1.1558  0. 0.  0.                           ! 18        
c 0.   0.    1.E6  -1.E6  1.E6  1.E6                                           
c 8. -1.                                       face 2                           
c 4  .1455   2.2670  -.6395  1.1558  0. 0.  0.                                   
c -18.6341335016   -18.6341335016   1.E6  -1.E6  1.E6  1.E6               ! 40   
c 0  0.      0.      0.      0.      0. 0.  0.                                   
c 0.  0.  1.E6  -1.E6  1.E6  1.E6 0.                                     ! 57             
c 2  64.                                                 IRD(=2, 25 or 4), resol(
c .1                                         step                               
c 2  315.5203145  -0.0698131700798  315.5203145  0.0698131700798        ! 63-67
c 'CHANGREF'
c XS -22.0023076169     ! RM*tan(AT-ACN+omga- == ACN-omga+ ==ACN)
c 'TRAROT'   MXQ1S01               SBEN             
c 0. 0. 0.   1.570796E+00 0. 0.

        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
      write(lw,fmt=txtfrm) '''DIPOLE''',name,' '//lbl2,it
      it = it + 1
c      at = abs(ang)*180.d0/4.d0*atan(1.d0) 
        write(lw,fmt='(I1,A)') i0,'  .Dipole'
        write(lw,fmt='(F14.7,2x,F15.8,T65,A)')  AT, RM,' ! AT, RM'
        write(lw,fmt='(F14.7,2x,2F15.8,a,T65,a)')acn,hnorm,b2,' 0. 0.'
     >  ,'        ! ACNT,  HNORM, indices'
        write(lw,fmt='(a)') ' 8. -1.         face 1 ' 
        write(lw,fmt='(a)') ' 4 .1455 2.2670 -.6395 1.1558 0. 0. 0.'
        write(lw,fmt='(2(F14.7,2x),a)')omgp,e1,' 1.E6 -1.E6 1.E6 1.E6'
        write(lw,fmt='(a)') ' 8. -1.         face 2 ' 
        write(lw,fmt='(a)') ' 4 .1455 2.2670 -.6395 1.1558 0. 0. 0.'
        e2deg = e2/(4.*atan(1.d0))*180.d0
        write(lw,fmt='(2(F14.7,2x),a)')omgm,-e2deg,
     >  ' 1.E6 -1.E6 1.E6 1.E6'
        write(lw,fmt='(a)') ' 0.  0.         face 3 ' 
        write(lw,fmt='(a)') ' 4 .1455 2.2670 -.6395 1.1558 0. 0. 0.'
        write(lw,fmt='(a)') ' 0. 0.   1.E6 -1.E6 1.E6 1.E6  0.'

        write(lw,fmt='(2A)') ' 2   64           '
        write(lw,fmt='(2A)') ' .2               '
        ffrad = ffang/180.*(4.*atan(1.d0))
        write(lw,fmt='(a,4(f13.7,1x),t70,A)') ' 2  '
     >  ,rers, -ffrad, rers, +ffrad, '  ! 63-67'

      write(lw,fmt=txtfrm) '''CHANGREF''',name,'      ',it
      it = it + 1
      write(lw,fmt='(a,f14.7,t65,a)') ' XS ',xs,
     >' ! RM*tan(ACN-omga+ ==ACN)'   

c         write(*,*) ' e2z i, txt4, var(i) // ',i, txt4,hnorm,b2,ReRs
c                   read(*,*)

      if(b<0.d0) 
     >write(lw,fmt=txtfrm) '''YMY''','YMY_OUT',' ',it

      if(Vbnd) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm) '''TRAROT''',name,' '//lbl2,it
        write(lw,fmt='(a,1p,e16.8,a)') '0. 0. 0. ',-tilt,' 0. 0.'
      endif

      goto 99

 4    continue
C----- QUAD
      if(xl .eq. 0.d0) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm) '''MARKER''',name,' '//lbl2,it
        goto 99
      endif
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,' '//lbl2,it
      write(lw,fmt='(I1,A)') i0,'  .Quad'
      write(lw,fmt='(F12.6,1x,F7.2,2(1x,F15.9),8(1x,F3.1))')
     >xl*cm,x10,x0,ak1*bro,x0,x0,x0,x0,x0,x0,x0,x0

      txt = txffq
      txt=txt//name
      if(frf .eq. 0.) txt = '0. 0. '//txt
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfq
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfq
      write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      write(lw,fmt='(2A)')  ' 0.5   ! cm ', kley
c      if(xl*cm-int(xl*cm).gt.0.5) then 
c            istepdip=int(xl*cm)+1
c      else
c            istepdip=int(xl*cm)
c      endif
c      if(istepdip .le. 0) istepdip = 2
c      write(lw,fmt='(A,I0,A,2X,A)') ' #30|',istepdip,'|30   Quad',name
c      write(lw,fmt='(A,2X,A)') '#10|10|10  Quad',name
      write(lw,fmt='(A)') '1 0. 0. 0.'
      goto 99

 5    continue
C----- SEXT ***
      if(xl .eq. 0.d0) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm) '''MARKER''',name,' '//lbl2,it
        goto 99
      endif
        it = it + 1
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,' '//lbl2,it
      write(lw,fmt='(I1,A)') i0,'  .Sext'
      b3 = 10.d0 * ak2 * bro *(x10/100.)**2 / 2.d0 
      write(lw,fmt='(F12.6,1x,F7.2,2(1x,F15.9),8(1x,F3.1))')
     >xl*cm,x10,x0,x0,b3,x0,x0,x0,x0,x0,x0,x0
      write(lw,fmt='(A)') '0.00 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') '0.00 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'

      write(lw,fmt='(2A)')  ' 0.5   ! cm ', kley

c      if(xl*cm-int(xl*cm).gt.0.5) then 
c            istepdip=int(xl*cm)+1
c      else
c            istepdip=int(xl*cm)
c      endif
c      if(istepdip.lt.10) then
c        write(lw,fmt='(2A)') ' #30|9|30   Sext',name
c      elseif(istepdip.ge.10.and.istepdip.lt.100) then
c        write(lw,fmt='(A,I2,A,2X,A)') ' #30|',istepdip,'|30   Sext',name
c      elseif(istepdip.ge.100.and.istepdip.lt.1000) then
c        write(lw,fmt='(A,I3,A,2X,A)') ' #30|',istepdip,'|30   Sext',name
c      elseif(istepdip.ge.1000) then
c        write(lw,fmt='(2A)') ' #30|9999|30   Sext',name
c      endif
C      write(lw,fmt='(A,2X,A)') '#10|10|10  Sext',name
      write(lw,fmt='(A)') '1 0. 0. 0.'
      goto 99

 6    continue
C----- OCTU ***
      if(xl .eq. 0.d0) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm) '''MARKER''',name,' '//lbl2,it
C        it = it - 1
        goto 99
      endif
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,' '//lbl2,it
      write(lw,fmt='(I1,A)') i0,'  .Octu'
      write(lw,fmt='(F12.6,1x,F7.2,2(1x,F15.9),8(1x,F3.1))')
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
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm)'''MARKER''',name,' '//lbl2,it
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
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,' '//lbl2,it
      write(lw,fmt='(I1,A)') i0,'  .Mult'
      write(lw,fmt='(F12.6,1x,F7.2,2(1x,F15.9),8(1x,F3.1))')
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
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
      write(lw,fmt=txtfrm) '''TRAROT''',name,' '//lbl2,it
      goto 99

 9    continue
C----- YROT ***
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
      write(lw,fmt=txtfrm) '''TRAROT''',name,' '//lbl2,it
      goto 99

 10   continue
C----- MARK, MATC(H)
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
      write(lw,fmt=txtfrm) '''MARKER''',name,' '//lbl2,it
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
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,' '//lbl2,it
      write(lw,fmt='(I1,A)') i0,'  .kicker'
      write(lw,fmt='(F12.6,1x,F7.2,2(1x,F15.9),7(1x,F3.1))')
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
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
      write(lw,fmt=txtfrm) '''DRIFT''',name,' '//lbl2,it
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
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm) '''MARKER''',name,' '//lbl2,it
C        it = it - 1
        goto 99
      endif
      r0 = 111.11
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
      write(lw,fmt=txtfrm) '''SOLENOID''',name,' '//lbl2,it
      write(lw,fmt='(I1,A)') i0,'  .soleno'
      write(lw,fmt='(1P,G14.7,G10.4,G14.7)') xl*cm,r0,e1*bro*t2kg
      write(lw,fmt='(2F10.2)') r0, r0
      write(lw,fmt='(A,2X,A)') '1.  Soleno',name
      write(lw,fmt='(A)') '1 0. 0. 0. '
      goto 99

 18   continue
C----- MATR
      if(xl .eq. 0.d0) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm)'''MARKER''',name,' '//lbl2,it
        goto 99
      else
        stop ' Element MATR not translated' 
      endif
      goto 99

 19   continue
C----- RFCA
      if(xl .eq. 0.d0) then
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
        write(lw,fmt=txtfrm)'''MARKER''',name,' '//lbl2,it
        goto 99
      endif
      if( okAutoref) then
       if(name .eq. 'R121' .or. name .eq. 'R221') then
        write(lw,fmt='(a)') '''AUTOREF'''
        write(lw,fmt='(a)') ' 1     ! Traj #1 is H-reference'
        write(lw,fmt='(a)') '''AUTOREF'''
        write(lw,fmt='(a)') ' 1.5   ! Traj #1 is H- & V-reference'
        write(lw,fmt='(a)') '''FAISCEAU'''
       endif
      endif
        lbl2 = fnwSeg(debstr(fnwSeg):finstr(fnwSeg)-4)
      write(lw,fmt=txtfrm) '''CAVITE''',name,' '//lbl2,it
      write(lw,fmt='(a)') '10     CBeta'
      write(lw,fmt='(1p,2(e23.15,2x))') xl, freq
C      write(lw,fmt='(1p,2(e23.15,2x),a)') volt, (phas-90.d0)/180.d0*pi, 
           phas = 1.53588957859 
      write(lw,fmt='(1p,2(e23.15,2x),a)') volt, phas, 
     >                           ' -2     ! 0 : boost+drift'
C     >                           ' 1     ! 0 : boost+drift'
C     >                           ' 0     ! 0 : boost+drift'
      goto 99

 22   continue
C----- MARK, PATC(H) ***
      Ys = xoff*1.d2
      Zs = yoff*1.d2
      Xs = zoff*1.d2
      ZR = xptch *180.d0/pi
        write(lw,fmt=txtfrm) '''CHANGREF''',name,'      ',it
        write(lw,fmt='(5(a,f14.7),a)') 
     >  'XS  0.  ! ZR ',zr,' XS ',xs,' YS ',Ys,' ZS ',Zs,' ! ( zr=',
     >   zr *pi/180.d0 , ' rad)'
      goto 99


 99    continue
c       write(*,fmt='(A,T12,A)') kley,name
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
      FUNCTION IDLUNI(
     >                LN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IDLUNI

      LOGICAL OPN

      I = 20
 1    CONTINUE
        INQUIRE(UNIT=I,ERR=99,IOSTAT=IOS,OPENED=OPN)
        I = I+1
        IF(I .EQ. 97) GOTO 99
        IF(OPN) GOTO 1
C        IF(IOS .GT. 0) GOTO 1
      
      LN = I-1
      IDLUNI = .TRUE.
      RETURN

 99   CONTINUE
      LN = 0
      IDLUNI = .FALSE.
      RETURN
      END
      FUNCTION GTTEXT(NRES,LUNR,TXT,
     >                              TXTIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL GTTEXT
      CHARACTER(*) TXT
      CHARACTER(*) TXTIN
      LOGICAL STRCON
      INTEGER DEBSTR, FINSTR

      READ(LUNR,FMT='(A)',ERR=99,END=98)  TXTIN

      DOWHILE(.NOT.STRCON(TXTIN,TXT(DEBSTR(TXT):FINSTR(TXT)),
     >                                                       IS))
        READ(LUNR,FMT='(A)',ERR=99,END=98)  TXTIN
      ENDDO

      GTTEXT = .TRUE.
      GOTO 10

 99   CONTINUE
      IF(NRES.GT.0) THEN
        WRITE(NRES,*) ' Pgm zgoubi/gttext - ERR upon read.' 
        WRITE(NRES,*) ' Text was : ',TXT(DEBSTR(TXT):FINSTR(TXT))
      ENDIF
      GTTEXT = .FALSE.
      GOTO 10

 98   CONTINUE
      IF(NRES.GT.0) THEN
        WRITE(NRES,*) ' Pgm zgoubi/gttext - EOF upon read.' 
        WRITE(NRES,*) ' Text was : ',TXT(DEBSTR(TXT):FINSTR(TXT))
      ENDIF
      GTTEXT = .FALSE.
      GOTO 10

 10   CONTINUE
      RETURN
      END
      subroutine getRijInd(name,
     >                       r11,r12,r33,r34)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IDLUNI, ok, strcon
C      character(*) name, matIndiv
      character(*) name
      character(800) txt800
      logical first
      integer debstr, finstr
      save lune
      data first / .true. /

      if(first) then
        ok = idluni(
     >            lune)
            write(*,*) ' ******',name


C        open(unit=lune,file=matIndiv)
        first = .false.
      endif

 1    continue
c        read(lune,fmt='(a)',err=99,end=98) txt800
        if(strcon(txt800,name(debstr(name):finstr(name)), 
     >                                              is)) then

c          read(txt800,*) s,r11,r12,r13,r14,r15,r16,
c     >    r21,r22,r23,r24,r25,r26,
c     >    r31,r32,r33,r34,r35,r36,
c     >    r41,r42,r43,r44,r45,r46
          r11=0.d0 ; r12=0.d0 ; r33=0.d0 ; r34=0.d0 

       else
         goto 1
       endif

      return

 97   stop 'Pgm getRijInd. Error during read txt800.'

 99   stop 'Pgm getRijInd. Error during read matIndiv.'
      
 98   write(*,*) ' End of file matIndiv.'
      write(*,*) ' was looking for this text : '
     >//name(debstr(name):finstr(name))
        stop 

      end 
