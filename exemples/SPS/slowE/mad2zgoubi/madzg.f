C23456789012345678901234567890123456789012345678901234567890123456789012
      implicit double precision (a-h,o-z)
      character fnr*80, fnw*80, fnw2*81, txt*132
      character txt590*590,keyword*20,txt20*20,txt40*40
      parameter(lr=9, lw=10, lw2=12, lout=14)
      parameter (mxl=15000)
      dimension iq(mxl)
      character kley*4, name*20
      parameter(mxk=24)
      character kle(mxk)*4, ny*1
      character*35 warn
      logical TOMANY

      logical ok, strcon
      integer debstr, finstr
      parameter (pi = 4.d0*atan(1.d0),dpi = 2.d0 * pi)

      data kle/ 'DRIF', 'RBEN', 'SBEN', 'QUAD', 'SEXT', 'OCTU', 'MULT',
     >'SROT', 'YROT', 'MARK', 'KICK', 'HMON', 'VMON','HKIC', 'VKIC', 
     >'SOLE', 'RCOL', 'MATR', 'MONI','TKIC','INST', 'HMON', 'VMON',
     >'RFCA' /

C This results from slightly wrong bend angle from MADX file. The compensaton is 
C re-distributed to the arc bends by the translator
      data devTot / dpi /   
C      data devTot / 6.2831853071795862d0 /    !  store
C      data devTot / 6.2783842424899960d0 /    ! injection
      data ksnk / 1 /

C      write(*,*) ' name of the survey file to be translated:'
C      read(*,fmt='(A)') fnr
C      write(*,*) ' name of the zgoubi data file:'
C      read(*,fmt='(A)') fnw

      fnw='trad.out'
      fnw2='trad.out2'

      fnr='madzg.in'
      if(fnr.ne.'madzg.in') then
        write(*,*)
        write(*,*) 'fnr MUST be  ''madzg.in'' '
        write(*,*)
      endif
      fnw='trad.out'

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
      
      write(*,*) 
      write(*,*) ' Fringe field coeff :',
     >             ' lhc/recycler/musr/default (1/2/3/99)'
C 2    read(*,*,err=2) kff
         kff = 99
C      if( kff .ne. 1 .and. kff .ne. 2 .and. kff .ne. 3) kff = 99
      write(*,*) ' Fringe field coeff : option set to ',kff
      call lmntff(kff)

      write(*,*) 
      write(*,*) ' Fringe fields on/off (1/0) :'
 3    read(*,*,err=3) frf
C          frf  = 0
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

      stp = 1.d0
      write(*,*)
      write(*,*) ' Step size in bends (cm, default is 1cm) :'
      read(*,*,err=9,end=9) stp
      goto 91
 9    continue
      stp = 1.d0
 91   continue
      if(stp.le.0.d0 .or. stp .ge. 100.) stp=1.d0
      call lmnstp(stp)

      write(*,*)
c      write(*,*) ' Include snake  (no/SPINR/fieldMap : 0/1/2) :'
c      read(*,*,err=8) ksnk
 8    continue
      if( ksnk .lt.0 .or. ksnk .gt. 2  ) ksnk = 1
c      write(*,*) ' ksnk = ',ksnk
c     read(*,*)
       ksnk = 0
      call lmnksn(ksnk)

C Get momentum
 32   continue
        read(lr,fmt='(A)',end=85,err=99) txt590
        ok =       strcon(txt590,'@ PC',
     >                                  IS) 
      if(.not. ok) goto 32        
      read(txt590(finstr(txt590)-17:finstr(txt590)),*) pmom

c       write(*,*) pmom       

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
c        write(*,*) txt590
c             read(*,*)

        read(txt590,fmt='(24e19.0,2e19.0,a,2e19.0)',end=85,err=99)
     >  xl, ang, ak1,ak2,ak3,ak4, tilt, e1,e2, h1,h2,sCntr,
     >  alfx,betx, alfy,bety, xco,yco, Dx,Dxp, Dy,Dyp, xmu,ymu, 
     >  hkic,vkic,
     >  txt40, xp, yp

        if(ir.eq.0) then
          it = 0
          pmom =  1.d0 * 2.99792458d8 / 1.d9
          call objet(lw,it,pmom*1.d9,xco,yco,xp,yp,
     >                           bro)        
          bro = 1.d0 

          write(*,*) ' bro, p : ',bro, pmom
        endif

        read(txt40,*,end=85,err=99) txt20,keyword
        kley = keyword
        name = txt20

C            if(hkic.ne.0.d0) then
c             write(*,*) trim(kley)
c            if(trim(kley) .eq. "VKIC") then
c        write(*,*)
c     >  xl, ang, ak1,ak2,ak3,ak4, tilt, e1,e2, h1,h2,sCntr,
c     >  alfx,betx, alfy,bety, xco,yco, Dx,Dxp, Dy,Dyp, xmu,ymu, 
c     >  hkic,vkic    ,txt40
c        write(*,*)
c     >  hkic,vkic    ,txt40
c             endif
     

c        write(*,*) '-----------------------------------'
c        write(*,*) kley, '  ', name, xl, ang, ak1, ak2
c             read(*,*)

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
        write(*,*) 'Press RETURN to continue '
        read(*,*)
        goto 86 
C        stop
 
 87     CONTINUE
        it = it + 1
        call lmnt
     >    (lw,lout,bro,frf,ik,noel,kley,name,devTot,xl,ang,
     >                     ak1,ak2,ak3,ak4,  tilt,e1,e2,h1,h2,
     >                                      hkic, vkic,
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
      write(lw,fmt='(A)')   ' '
      write(lw,fmt='(A,T111,I6)') '''OPTIONS''',it
      write(lw,fmt='(A)')   ' 1 1 '
      write(lw,fmt='(A)')   ' WRITE ON '
      write(lw,fmt='(A)')   ' '

C TWISS w/ search orbit ---------
      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''FIT''',it
      write(lw,fmt='(A)')   ' 2'
      write(lw,fmt='(A)')   ' 1  30 0  [-1.,1.]         !  Y_0'
      write(lw,fmt='(A)')   ' 1  31 0  [-1.,1.]         !  T_0'
      write(lw,fmt='(A)')   ' 2   1e-10  9999'
      write(lw,fmt='(A)')   ' 3.1 1 2 #End .0   1. 0     ! Y_co'
      write(lw,fmt='(A)')   ' 3.1 1 3 #End .0   1. 0     ! T_co'
C-------------------------------
      write(lw,fmt='(A)') ' '
      write(lw,fmt='(A,T111,I6)') '! ''MATRIX''',it
      write(lw,fmt='(A)')   '! 1  11 '
      write(lw,fmt='(A)') ' '
      write(lw,fmt='(A,T111,I6)') '!  ''FIT''',it
      write(lw,fmt='(a)') '! 4  ' 
      write(lw,fmt='(a)') '! 1 30 0 [-1.,1.]    !   Y_0 '
      write(lw,fmt='(a)') '! 1 31 0 [-1.,1.]    !   T_0 '
      write(lw,fmt='(a)') '! 5  20 0 .1         !   QF '
      write(lw,fmt='(a)') '! 5  24 0 .1         !   QD '
      write(lw,fmt='(a)') '! 4   1e-10  9999'
      write(lw,fmt='(a)') '! 3.1 1 2 #End 0. 1. 0           ! Y_co'
      write(lw,fmt='(a)') '! 3.1 1 3 #End 0. 1. 0           ! T_co'
      write(lw,fmt='(a)') '! 0.1 7 7 #End 0.666666 1. 0     ! Q_Y'
      write(lw,fmt='(a)') '! 0.1 8 8 #End 0.58     1. 0     ! Q_Z'
      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''TWISS''',it
      write(lw,fmt='(A)')   '2  1. 1. '
      write(lw,fmt='(A)') ' '
      write(lw,fmt='(A)') ' '
      write(lw,fmt='(A)') '''END'''
      it = it + 1
      write(lw,fmt='(A)') ' '
      write(lw,fmt='(A)') '''REBELOTE'''
      write(lw,fmt='(A)') ' 49999 0.1 99'
      write(lw,fmt='(A,T111)') '''FAISCEAU'''
      write(lw,fmt='(A)') ' '
      write(lw,fmt='(A)') '''END'''

      write(*,*) ' Read ',ir,' elements from madzg.in file ',fnr
      write(*,*) ' Translated to',it,
     >' elements into the zgoubi data file ', fnw
      write(*,*) ' end of madzg.in file '

      call lmnAle(
     >            aleTot, nbArcBnd)

      write(*,*) ' '
      write(*,*) '-----------------------------------------------'
      write(*,*) '-----------------------------------------------'
      write(*,*) 'Total deviation in BENDs with KPOS=3, ALE_tot = '
     >,-2.d0*aleTot,'    (expected is 2pi = ',dpi,')'
      write(*,*) '(found ',nbArcBnd,' arc bends)'
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) 'If not 2*pi, recompile madzg with devTot=2pi '
     >//'in madzg.f, re-run, this will give the value '
     >//'(ALE_tot) to be injected in devTot.'
      write(*,*) '-----------------------------------------------'
      write(*,*) '-----------------------------------------------'
      write(*,*) ' '

      goto 999

 97   write(*,*) ' error open madzg.in file'
      goto 999
 99   write(*,*) ' error during read of madzg.in file'
      goto 999

 999  continue
      close(lr)
      close(lw)
      close(lw2)
      stop
      end
        
      subroutine objet(lw,it,pmom,xco,yco,xpco,ypco,
     >                            bro)
      implicit double precision (a-h,o-z)
     
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
      write(lw,fmt='(A)') '5' 
      write(lw,fmt='(A)') '.001 .001 .001 .001 0. .0001  '
      write(lw,fmt='(A)') '.0 .0 .0 .0 0.  1.  '
c      write(lw,fmt='(1p,4(e14.6,2x),A)') 
c     >xco*1d2,xpco*1d3,yco*1d2,ypco*1d3, ' 0. 1.'


      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''PARTICUL''',it
      write(lw,FMT='(A)') '9.3827203E+02 1.602176487E-19 1.7928474 0 0' 

      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''FAISCEAU''',it
      it = it + 1
      write(lw,fmt='(A,T111,I6)') '''FAISTORE''',it
      write(lw,fmt='(A)') 'zgoubi.fai    E-Septum_up'
      write(lw,fmt='(A)') '1'

c      it = it + 1
c      write(lw,fmt='(A,T111,I6)') '''OPTICS''',it
c      write(lw,FMT='(A)') 
c     > ' 2   Print out transport coeffs to zgoubi_MATRIX_out' 

      it = it + 1
      write(lw,FMT='(A)') ' ' 
      write(lw,fmt='(A,T111,I6)') '''SCALING''',it 
      write(lw,fmt='(A)') '1  12 '
      write(lw,fmt='(A)') 'BEND' 
      write(lw,fmt='(A)') '   -1   ' 
      write(lw,fmt='(1P,E16.8)') BRO
      write(lw,fmt='(A)') '1  '
      write(lw,fmt='(A)') 'MULTIPOL '
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8)') BRO
      write(lw,fmt='(A)') '1      '
      write(lw,fmt='(A)') 'MULTIPOL HKIC'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8)') BRO
      write(lw,fmt='(A)') '1      '
      write(lw,fmt='(A)') 'MULTIPOL VKIC'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8)') BRO
      write(lw,fmt='(A)') '1      '
      write(lw,fmt='(A)') 'MULTIPOL QF*'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8,T20,A)') BRO*1.0005153,' ! FIT# 20'
      write(lw,fmt='(A)') '1      '
      write(lw,fmt='(A)') 'MULTIPOL QD*'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8,T20,A)') BRO*0.99939388,' ! FIT# 24'
      write(lw,fmt='(A)') '1      '
      write(lw,fmt='(A)') 'MULTIPOL '
     >//'LSE.10602 LSE.22402 LSE.40602 LSE.52402'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1p,E16.8A)') 0.,'   ! * (-0.011992)'
      write(lw,fmt='(A)') '1      '
      write(lw,fmt='(A)') 'MULTIPOL LSFA* LSFC*'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8,T20,A)') BRO*1.8280002,' ! FIT# 32'
      write(lw,fmt='(A)') '1      '
      write(lw,fmt='(A)') 'MULTIPOL LSFB*'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8,T20,A)') BRO*0.32674538,' ! FIT# 36'
      write(lw,fmt='(A)') '1      '
      write(lw,fmt='(A)') 'MULTIPOL LSDA*'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8,T20,A)') BRO*0.94126777,' ! FIT# 40'
      write(lw,fmt='(A)') '1      '
      write(lw,fmt='(A)') 'MULTIPOL LSDB*'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8,T20,A)') BRO*1.0000004,' ! FIT# 44'
      write(lw,fmt='(A)') '1      '
      write(lw,fmt='(A)') 'MULTIPOL RB_*'
      write(lw,fmt='(A)') '-1         '
      write(lw,fmt='(1P,E16.8,A)')
     >             BRO* (+0.374),' ! 3.74cm at max, 3.61cm E-septum' 
      write(lw,fmt='(A)') '1      '

      it = it + 1
      write(lw,fmt='(A)')   ' '
      write(lw,fmt='(A,T111,I6)') '''OPTIONS''',it
      write(lw,fmt='(A)')   ' 1 1 '
      write(lw,fmt='(A)')   ' WRITE ON '
      write(lw,fmt='(A)')   ' '

      return
      end 

      subroutine lmnt
     >(lw,lout,bro,frf,ik,noel,kley,name,devTot,xl,ang,
     >           ak1,ak2,ak3,ak4,    tilt,e1,e2,h1,h2,
     >                                      hkic, vkic,
     >                                it)
      implicit double precision (a-h,o-z)
      character kley*(*)

      character(20) name
      character ny*1

      character*80 txfd, txfq,fqlhc , fdlhc ,fqrec , fdrec,fd,fq
      character*80 txffd,txffq,ffqlhc,ffdlhc,ffqrec,ffdrec,ffd,ffq
      character*80 ffdmu, fdmu, ffqmua,ffqmus,fqmua
      save txfd,txfq,txffd, txffq
      logical drimul
      save drimul
      logical strcon, ok

      parameter (cm=100.d0, t2kg=10.d0)
      parameter(i0=0, i1=1, i2=2, i4=4, i6=6)
      parameter(x0=0.d0,x1=1.d0,x2=2.d0,x3=3.d0,x4=4.d0,x5=5.d0,x6=6.d0)
      parameter(x7=7.d0,x8=8.d0,x10=10.d0,x20=20.d0,
     >                                    x30=30.d0,x999=999.d0)

      character txt*80, txtfrm*23
      save kpos, kbm
      save aleTot, nbArcBnd, ksnk

      save stpBnd

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
     
              data numb / 0 /
       data stpBnd / 1.d0 /

c      parameter (stpsiz = 2.d0)  ! step size in cm
      parameter (stpsiz = 1.d0)  ! step size in cm
c      parameter (stpsiz = 0.5d0)  ! step size in cm
c      parameter (stpsiz = 0.25d0)  ! step size in cm

      pi = 4.d0*atan(1.d0)
      zero=0.d0
      txtfrm = '(A,T12,A,T37,A,T111,i6)'

C----- 1    2    3    4    5    6    7    8    9    10   11   12   13
C----- DRIF RBEN SBEN QUAD SEXT OCTU MULT SROT YROT MARK KICK HMON VMON
C----- 14   15   16   17   18   19   20   21   22   23   24
C----- HKIC VKIC SOLE RCOL MATR MONI TKIC INST HMON VMON RFCA
       goto (
     > 1,   2 ,  2 ,  4,   5,   6,   7,   8,   9,   10,  11,  12,  12,
     > 11,  11,  13,  1,  18,   12, 1, 1, 12, 12, 1 ) ik

 1    continue
C----- DRIF      
      if(xl .eq. 0.d0) then
        write(lw,fmt=txtfrm)'''MARKER''',name,kley,it
        goto 99
      else
          if(drimul) then
            dum=1.d-20
            write(lw,fmt=txtfrm) '''MULTIPOL''',kley,name,it
            write(lw,fmt='(I1,A)') i0,'  .Drift'
            write(lw,fmt='(G12.6,1X,G8.2,1X,G16.10,1X,9(G7.1,1X))')
     >      xl*cm,x10,x0,dum,x0,x0,x0,x0,x0,x0,x0,x0
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
C case RBEN (madx twiss writes xl as arc length)
        xlarc = xl   !!  /(2.d0*sin(ang/2.d0)) * ang
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
        fac = 1.d-20
      endif
      write(lw,fmt='(F11.6,1X,F7.2,1X,3(F13.8,1X),7(F5.3,1X))')
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
      if(xlmag*cm/stpsiz-int(xlmag*cm/stpsiz) .gt. 0.5d0) then 
            istepdip=int(xlmag*cm/stpsiz)+1
      else
            istepdip=int(xlmag*cm/stpsiz)
      endif
      if(istepdip.lt.10) then
        write(lw,fmt='(2A)') ' #30|9|30    Dip',name
      elseif(istepdip.ge.10.and.istepdip.lt.100) then
        write(lw,fmt='(A,I2,A,2X,A)') ' #30|',istepdip,'|30    Dip',name
      elseif(istepdip.ge.100.and.istepdip.lt.1000) then
        write(lw,fmt='(A,I3,A,2X,A)') ' #30|',istepdip,'|30    Dip',name
      elseif(istepdip.ge.1000.and.istepdip.lt.10000) then
        write(lw,fmt='(A,I4,A,2X,A)') ' #30|',istepdip,'|30    Dip',name
      elseif(istepdip.ge.10000) then
        write(lw,fmt='(2A)') ' #30|99999|30    Dip',name
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
C-----RBEN -> BEND
        if( 
     >  trim(name) .eq. 'MPSH_RB.21202' .or. 
     >  trim(name) .eq. 'MPLH_RB.21431' .or. 
     >  trim(name) .eq. 'MPNH_RB.21732' .or. 
     >  trim(name) .eq. 'MPLH_RB.21995' .or. 
     >  trim(name) .eq. 'MPLH_RB.22195' 
     >                         ) then
          if( 
     >    trim(name) .eq. 'MPSH_RB.21202' ) then
            b1 = 7.6765e-5
          elseif( 
     >    trim(name) .eq. 'MPLH_RB.21431' ) then
            b1 = -0.49e-3
          elseif(  
     >    trim(name) .eq. 'MPNH_RB.21732' ) then
            b1 = -0.33309e-3
          elseif( 
     >    trim(name) .eq. 'MPLH_RB.21995' ) then
            b1 = -0.2503e-3
          elseif( 
     >    trim(name) .eq. 'MPLH_RB.22195'  ) then
            b1 = .35585e-3
          endif
          
           write(lw,fmt=txtfrm) '''MULTIPOL''',
     >     'RB_'//name(1:4)//name(8:13),'bump_RB',it
          write(lw,fmt='(I1,A)') i0,'  bump_RB'
          write(lw,fmt='(F11.6,1X,F6.2,1X,F16.10,A)')
     >    xl*cm, x10,b1*bro*T2kG,' 0. 0. 0. 0. 0. 0. 0. 0. 0. '
          txt=txffd
          txt='0. 0. '//txt
          write(lw,fmt='(A)') txt
          write(lw,fmt='(A)') txfd
          write(lw,fmt='(A)') txt
          write(lw,fmt='(A)') txfd
          write(lw,fmt='(1p,A)') 
     >    '0. 0. 0. 0.  0. 0. 0. 0. 0. 0.'
          write(lw,fmt='(2A)') ' #30|11|30    bump_RB',name
          write(lw,fmt='(A)') ' 2 0. 0. 0. '

          goto 99
          
        endif
      
      if(xl * ang .eq. 0.d0) then
        if(xl .eq. 0.d0) then
          write(lw,fmt=txtfrm) '''MARKER''',kley,name,it
        else
          write(lw,fmt=txtfrm) '''DRIFT''',kley,name,it
          write(lw,fmt='(1P,E14.6)') xl*cm
        endif
        goto 99
      endif

        ok = strcon(name,'_DH',
     >                         IS)
        ok = ok .and. name(1:2) .eq. 'BO'
     >          .or.  name(1:2) .eq. 'BI'
        if(ok) then
          read(name(is+3:is+4),*) num
          if(num .ge. 10 .and. num .le. 21) then

               numb = numb + 1
            write(86,*)  ' numb ',ang, (2.d0*pi - devTot)/dble(132),numb

            ang = ang + (2.d0*pi - devTot)/dble(132)
          
          endif
        endif

      write(lw,fmt=txtfrm) '''BEND''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .Bend'
      ro=xl/2./sin(ang/2.)
      b=ang * bro / xl *10.  !bro/ro*t2kg
C      xxl=ro*ang*cm
      ro=ro*cm

C err. FM June 2009      write(lw,fmt='(2F14.7,F15.8)') xl,tilt,b 
      if(abs(tilt).lt.1d-10) then 
        fac = 1.d0
      else
        fac = 1.d-20
      endif
      write(lw,fmt='(2F16.9,F15.8)') xl*cm,tilt,fac*b 
      te=ang/2.+e1
      if(frf .eq. 0.) then
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,te
      else 
C        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,te
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x30,x10,te
      endif
      write(lw,fmt='(A)') txfd
      ts=ang/2.+e2
      if(frf .eq. 0.) then
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,ts
      else
C        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,ts
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x30,x10,te
      endif
      write(lw,fmt='(A)') txfd
C      step = xl
C      step = stpsiz
      step = stpBnd      
      write(lw,fmt='(1P,E12.4,A)') dble(nint(step)),'  Bend'
C      write(lw,fmt='(F12.6,A4)') step * 10.,'  Bend'
C      write(lw,fmt='(A,2X,A)') ' #10|10|10   Bend',name
      if(kpos.eq.3) then
        if(abs(tilt).lt.1d-10) then 
          write(lw,fmt='(A)') '3 0. 0. 0.'
        else
          write(lw,fmt='(A)') '1 0. 0. 0.'
        endif
        aleTot = aleTot - ang/2.d0
        ok = strcon(name,'_DH',
     >                         IS)
        ok = ok .and. name(1:2) .eq. 'BO'
     >          .or.  name(1:2) .eq. 'BI'
        if(ok) then
c          write(*,*) ' name   ***',name,'***'
c          write(*,*) ' name   ***',name(is+3:is+4),'***'
c          read(*,*)
          read(name(is+3:is+4),*) num
          if(num .ge. 10 .and. num .le. 21) then
            nbArcBnd = nbArcBnd + 1
          endif
        endif
      else
        write(lw,fmt='(1P,I2,3(1x,e18.10))') kpos,xce,yce,ale
      endif
      goto 99

 3    continue
C----- SBEN -> BEND
      if(xl * ang .eq. 0.d0) then
        if(xl .eq. 0.d0) then
          write(lw,fmt=txtfrm) '''MARKER''',kley,name,it
        else
          write(lw,fmt=txtfrm) '''DRIFT''',kley,name,it
          write(lw,fmt='(1P,E14.6)') xl*cm
        endif
        goto 99
      endif

        ok = strcon(name,'_DH',
     >                         IS)
        ok = ok .and. name(1:2) .eq. 'BO'
     >          .or.  name(1:2) .eq. 'BI'
        if(ok) then
          read(name(is+3:is+4),*) num
          if(num .ge. 10 .and. num .le. 21) then

               numb = numb + 1
            write(86,*)  ' numb ',ang, (2.d0*pi - devTot)/dble(132),numb

            ang = ang + (2.d0*pi - devTot)/dble(132)
          endif
        endif

      write(lw,fmt=txtfrm) '''BEND''',kley,name,it
      write(lw,fmt='(I1,A)') i0,'  .Bend'
      ro=xl/ang
      xxl = 2.d0*ro*sin(ang/2.d0)
      b=bro/ro*t2kg
              write(lout,*) ' ********* tilt ',tilt
      if(abs(tilt).lt.1d-10) then 
        fac = 1.d0
      else
        fac = 1.d-20
      endif
      write(lw,fmt='(1P,3(1X,E17.9))') xxl*cm,tilt,fac*b 
      te=e1
      if(frf .eq. 0.) then
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,te
      else
C        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,te
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x30,x10,te
      endif
      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
      ts=e2
      if(frf .eq. 0.) then
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x0,x0,ts
      else
C        write(lw,fmt='(2F6.2,2(1x,F12.8))') x20,x8,ts
        write(lw,fmt='(2F6.2,2(1x,F12.8))') x30,x10,te
      endif
      write(lw,fmt='(A)') '4 .2401  1.8639  -.5572  .3904 0. 0. 0.'
C      istepdip=int(xxl*cm/3.0d0)
C      istepdip=int(xxl*cm/stpsiz)
      istepdip=int(xxl*cm/stpBnd)
      if(istepdip .ge. 10000) then
        if(istepdip .ge. 99999) istepdip = 99999
        write(lw,fmt='(A,I5,A,2X,A)') '#200|',istepdip,
     + '|200    Bend',name
      elseif(istepdip .ge. 1000  .and. istepdip .lt. 10000) then
        write(lw,fmt='(A,I4,A,2X,A)') '#200|',istepdip,
     + '|200    Bend',name
      elseif(istepdip .ge. 100  .and. istepdip .lt. 1000) then
        write(lw,fmt='(A,I3,A,2X,A)') '#200|',istepdip,
     + '|200    Bend',name
      elseif(istepdip .ge. 10  .and. istepdip .lt. 100) then
        write(lw,fmt='(A,I2,A,2X,A)') '#200|',istepdip,
     + '|200    Bend',name
      elseif(istepdip .lt. 10) then
        write(lw,fmt='(2A)') '#200|5|200    Bend',name 
      endif
      if(kpos.eq.3) then
        if(abs(tilt).lt.1d-10) then 
          write(lw,fmt='(A,1p,2x,e14.6)') '3 0. 0. ',-ang/2.d0
        else
          write(lw,fmt='(A)') '1 0. 0. 0.'
        endif
        aleTot = aleTot - ang/2.d0
        ok = strcon(name,'_DH',
     >                         IS)
        ok = ok .and. name(1:2) .eq. 'BO'
     >          .or.  name(1:2) .eq. 'BI'
        if(ok) then
c          write(*,*) ' name   ***',name,'***'
c          write(*,*) ' name   ***',name(is+3:is+4),'***'
c          read(*,*)
          read(name(is+3:is+4),*) num
          if(num .ge. 10 .and. num .le. 21) then
            nbArcBnd = nbArcBnd + 1
          endif
        endif
      else
        write(lw,fmt='(1P,I2,3(1x,e18.10))') kpos,xce,yce,ale
      endif
      goto 99

 4    continue
C----- QUAD
      if(xl .eq. 0.d0) then
        write(lw,fmt=txtfrm) '''MARKER''',kley,name,it
        goto 99
      endif
      call quanam(
     >             name)
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .Quad'
      write(lw,fmt='(F12.6,1X,F6.2,1X,2(F16.10,1X),8(F4.1,1X))')
     >xl*cm,x10,x0,ak1/xl*bro,x0,x0,x0,x0,x0,x0,x0,x0
      txt = txffq
      txt=txt//name
!      if(frf .eq. 0.) txt = '0. 0. '//txt
      txt = '0. 0. '//txt
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfq
      write(lw,fmt='(A)') txt
      write(lw,fmt='(A)') txfq
      write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      if(xl*cm/stpsiz-int(xl*cm/stpsiz).gt.0.5) then 
            istepdip=int(xl*cm/stpsiz)+1
      else
            istepdip=int(xl*cm/stpsiz)
      endif
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
      call sxtnam(
     >          name)
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .Sext'
      b3 = 10.d0 * ak2/xl * bro *(x10/100.)**2 / 2.d0
      if(
     >name .eq. 'LSE.10602'  .or.
     >name .eq. 'LSE.22402'  .or.
     >name .eq. 'LSE.40602'  .or.
     >name .eq. 'LSE.52402'  ) then
         b3 = 1.d0
      endif
      write(lout,*) ' SEXT         b3 = ',b3,' kG'
      write(lw,fmt=
     >'(F12.6,1X,F6.2,1X,2(F10.4,1X),F16.10,1X,7(F4.1,1X))')
     >xl*cm,x10,x0,x0,b3,x0,x0,x0,x0,x0,x0,x0
      write(lw,fmt='(A)') '0.00 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') '0.00 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      if(xl*cm/stpsiz-int(xl*cm/stpsiz).gt.0.5) then 
            istepdip=int(xl*cm/stpsiz)+1
      else
            istepdip=int(xl*cm/stpsiz)
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
      write(lw,fmt='(F10.4,1x,F6.2,1x,1p,4(E12.4,1x),6(E8.2,1x))')
     >xl*cm,x10,x0,x0,x0,x0,x0,x0,x0,x0,x0,x0
      write(lw,fmt='(A)') '0.000 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
      write(lw,fmt='(A)') '0.000 0.00  1.00 0.0 .0 0.0 .0 0.0 .0 0.0 .0'
      write(lw,fmt='(A)') txfq//name
        write(lw,fmt='(A)') ' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.'
      if(xl*cm/stpsiz-int(xl*cm/stpsiz).gt.0.5) then 
            istepdip=int(xl*cm/stpsiz)+1
      else
            istepdip=int(xl*cm/stpsiz)
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
      write(lw,fmt=txtfrm) '''MULTIPOL''',name,kley,it
      write(lw,fmt='(I1,A)') i0,'  .Mult'
      write(lw,fmt='(F9.4,1x,F6.2,1x,1p,3(e11.4,1x),7(e8.1,1x))')
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
      if(xl*cm/stpsiz-int(xl*cm/stpsiz).gt.0.5) then 
            istepdip=int(xl*cm/stpsiz)+1
      else
            istepdip=int(xl*cm/stpsiz)
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
      if    (trim(name) .eq. 'AP.UP.ZS21633') then
        write(lw,fmt='(a)') '''COLLIMA''  E-Septum_up'
        write(lw,fmt='(a)') ' 2'
        write(lw,fmt='(a)') ' 1.1  -20. 6.8  -20. 20.'
      endif
      it = it - 1
      write(lw,fmt=txtfrm) '''MARKER''',name,kley,it


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
c             write(*,*) ' hkic, vkic :', hkic, vkic,bro,b,t2kg
c                read(*,*)
      write(lw,fmt=txtfrm) '''MULTIPOL''',kley,name,it
      write(lw,fmt='(I1,A)') i0,'  .kicker'
      write(lw,fmt='(E14.6,1x,F6.2,1x,E14.6,1x,F12.6,1x,8(F4.1,1x))')
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
      write(lw,fmt=txtfrm) '''DRIFT''',name,kley,it
      write(lw,fmt='(F12.4)') xl*cm
c      write(*,*) '**',trim(name),'**',trim(name) .eq. 'BI9_B7.1',ksnk
      if    (trim(name) .eq. 'BI9_B7.1') then
        it = it+1
        write(lw,fmt='(A,T12,A,T111,i6)') '''MARKER''','SNAKE1',it
        if    (ksnk .eq. 1) then
          it = it+1
          write(lw,fmt='(A,T12,A,T111,i6)') '''SPINR''','SNK1',it
          write(lw,fmt=txtfrm) ' 1'
          write(lw,fmt=txtfrm) ' -45.   180.    snake_1'
        elseif(ksnk .eq. 2) then
              stop ' ksnk =2 NOT IMPLEMENTED'
          it = it+1
          write(lw,fmt=txtfrm) '''DRIFT''',' DRIF ','DRIFT_81',it 
          write(lw,fmt='(a)')   '  -655.441800'
c------------- INSERT SNAKE MAPS HERE
C          it = it+1
c          write(lw,fmt=txtfrm) '''DRIFT''',' VMON ','BO3_B7.1',it 
c          write(lw,fmt=txtfrm)        0.0000
C------------------
          it = it+1
          write(lw,fmt=txtfrm) '''DRIFT''',' DRIF ','DRIFT_81',it 
          write(lw,fmt='(a)') ' -655.441800'
        endif
      elseif(trim(name) .eq. 'BO3_B7.1') then
        it = it+1
        write(lw,fmt='(A,T12,A,T111,i6)') '''MARKER''','SNAKE2',it
        if    (ksnk .eq. 1) then
          it = it+1
          write(lw,fmt='(A,T12,A,T111,i6)') '''SPINR''','SNK2',it
          write(lw,fmt=txtfrm) ' 1'
          write(lw,fmt=txtfrm) ' +45.   180.    snake_2'
        elseif(ksnk .eq. 2) then
        endif
      endif
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
      r0 = 111.11
      write(lw,fmt=txtfrm) '''SOLENOID''',kley,name,it
      write(lw,fmt='(I1,A)') i0,'  .soleno'
      write(lw,fmt='(1P,G14.7,G10.4,G14.7)') xl*cm,r0,e1*bro*t2kg
      write(lw,fmt='(2F10.2)') r0, r0
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
c       write(6,fmt='(A,T12,A)') kley,name
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
      aleTot = 0.d0
      nbArcBnd = 0 
      return

      entry lmnAle(
     >             aleTotO,nbArcBndO)
      aleTotO = aleTot 
      nbArcBndO =nbArcBnd
      return

      entry lmndri(ny)
      drimul=ny.eq.'y'
      return

      entry lmnkbm(kbmi)
      kbm =kbmi 
      return

      entry lmnstp(stpi)
      stpBnd = stpi
      return

      entry lmnksn(ksnki)
      ksnk = ksnki
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
      SUBROUTINE SXTNAM(
     >                name)
      character(*) name
      if    (name(1:3) .eq. 'LSF') then
        if(name(5:9) .eq. '10805' .or.
     >      name(5:9) .eq. '11205' .or.
     >      name(5:9) .eq. '12005' .or.
     >      name(5:9) .eq. '12405' .or.
     >      name(5:9) .eq. '13205' .or.
     >      name(5:9) .eq. '13605' .or.
     >      name(5:9) .eq. '20805' .or.
     >      name(5:9) .eq. '21205' .or.
     >      name(5:9) .eq. '22005' .or.
     >      name(5:9) .eq. '22405' .or.
     >      name(5:9) .eq. '23205' .or.
     >      name(5:9) .eq. '23605' .or.
     >      name(5:9) .eq. '30805' .or.
     >      name(5:9) .eq. '31205' .or.
     >      name(5:9) .eq. '32005' .or.
     >      name(5:9) .eq. '32405' .or.
     >      name(5:9) .eq. '43205' .or.
     >      name(5:9) .eq. '53205' .or.
     >      name(5:9) .eq. '60805' .or.
     >      name(5:9) .eq. '61205' .or.
     >      name(5:9) .eq. '62005' .or.
     >      name(5:9) .eq. '62405' .or.
     >      name(5:9) .eq. '63205' .or.
     >      name(5:9) .eq. '63605' ) then
             name= 'LSFA'//name(4:20)
        elseif(   name(5:9) .eq. '10205' .or.
     >      name(5:9) .eq. '11405' .or.
     >      name(5:9) .eq. '12605' .or.
     >      name(5:9) .eq. '20205' .or.
     >      name(5:9) .eq. '21405' .or.
     >      name(5:9) .eq. '22605' .or.
     >      name(5:9) .eq. '30205' .or.
     >      name(5:9) .eq. '31405' .or.
     >      name(5:9) .eq. '32605' .or.
     >      name(5:9) .eq. '40205' .or.
     >      name(5:9) .eq. '41405' .or.
     >      name(5:9) .eq. '42605' .or.
     >      name(5:9) .eq. '50205' .or.
     >      name(5:9) .eq. '51405' .or.
     >      name(5:9) .eq. '52605' .or.
     >      name(5:9) .eq. '60205' .or.
     >      name(5:9) .eq. '61405' .or.
     >      name(5:9) .eq. '62605' ) then
             name= 'LSFB'//name(4:20)
        elseif(   name(5:9) .eq. '33205' .or.
     >      name(5:9) .eq. '33605' .or.
     >      name(5:9) .eq. '40805' .or.
     >      name(5:9) .eq. '41205' .or.
     >      name(5:9) .eq. '42005' .or.
     >      name(5:9) .eq. '42405' .or.
     >      name(5:9) .eq. '43605' .or.
     >      name(5:9) .eq. '50805' .or.
     >      name(5:9) .eq. '51205' .or.
     >      name(5:9) .eq. '52005' .or.
     >      name(5:9) .eq. '52405' .or.
     >      name(5:9) .eq. '53605' .or.
     >      name(5:9) .eq. '60805' .or.
     >      name(5:9) .eq. '61205' ) then
             name= 'LSFC'//name(4:20)
         endif
      elseif(name(1:3) .eq. 'LSD') then
        if(name(5:9) .eq. '10105' .or.
     >   name(5:9) .eq. '11305' .or.
     >   name(5:9) .eq. '12505' .or.
     >   name(5:9) .eq. '20105' .or.
     >   name(5:9) .eq. '21305' .or.
     >   name(5:9) .eq. '22505' .or.
     >   name(5:9) .eq. '30105' .or.
     >   name(5:9) .eq. '31305' .or.
     >   name(5:9) .eq. '32505' .or.
     >   name(5:9) .eq. '40105' .or.
     >   name(5:9) .eq. '41305' .or.
     >   name(5:9) .eq. '42505' .or.
     >   name(5:9) .eq. '50105' .or.
     >   name(5:9) .eq. '51305' .or.
     >   name(5:9) .eq. '52505' .or.
     >   name(5:9) .eq. '60105' .or.
     >   name(5:9) .eq. '61305' .or.
     >   name(5:9) .eq. '62505' ) then
             name= 'LSDA'//name(4:20)        
        elseif(name(5:9) .eq. '10305' .or.
     >   name(5:9) .eq. '11105' .or.
     >   name(5:9) .eq. '11505' .or.
     >   name(5:9) .eq. '12305' .or.
     >   name(5:9) .eq. '12705' .or.
     >   name(5:9) .eq. '13505' .or.
     >   name(5:9) .eq. '20305' .or.
     >   name(5:9) .eq. '21105' .or.
     >   name(5:9) .eq. '21505' .or.
     >   name(5:9) .eq. '22305' .or.
     >   name(5:9) .eq. '22705' .or.
     >   name(5:9) .eq. '23505' .or.
     >   name(5:9) .eq. '30305' .or.
     >   name(5:9) .eq. '31105' .or.
     >   name(5:9) .eq. '31505' .or.
     >   name(5:9) .eq. '32305' .or.
     >   name(5:9) .eq. '32705' .or.
     >   name(5:9) .eq. '33505' .or.
     >   name(5:9) .eq. '40305' .or.
     >   name(5:9) .eq. '41105' .or.
     >   name(5:9) .eq. '41505' .or.
     >   name(5:9) .eq. '42305' .or.
     >   name(5:9) .eq. '42705' .or.
     >   name(5:9) .eq. '43505' .or.
     >   name(5:9) .eq. '50305' .or.
     >   name(5:9) .eq. '51105' .or.
     >   name(5:9) .eq. '51505' .or.
     >   name(5:9) .eq. '52305' .or.
     >   name(5:9) .eq. '52705' .or.
     >   name(5:9) .eq. '53505' .or.
     >   name(5:9) .eq. '60305' .or.
     >   name(5:9) .eq. '61105' .or.
     >   name(5:9) .eq. '61505' .or.
     >   name(5:9) .eq. '62305' .or.
     >   name(5:9) .eq. '62705' .or.
     >   name(5:9) .eq. '63505'  ) then
             name= 'LSDB'//name(4:20)                  
        endif   
      endif
      return
      end
      SUBROUTINE QUANAM(
     >                name)
      character(*) name
      if    (name(1:3) .eq. 'QF.') then
        if(
     >    name(4:8) .eq. '31810' .or. 
     >    name(4:8) .eq. '32010' .or. 
     >    name(4:8) .eq. '32210' .or. 
     >    name(4:8) .eq. '32410' .or. 
     >    name(4:8) .eq. '32610' .or. 
     >    name(4:8) .eq. '32810' .or. 
     >    name(4:8) .eq. '33010' .or. 
     >    name(4:8) .eq. '33210' .or. 
     >    name(4:8) .eq. '33410' .or. 
     >    name(4:8) .eq. '40010' .or. 
     >    name(4:8) .eq. '40210' .or. 
     >    name(4:8) .eq. '40410' .or. 
     >    name(4:8) .eq. '40610' .or. 
     >    name(4:8) .eq. '40810' .or. 
     >    name(4:8) .eq. '41010' .or. 
     >    name(4:8) .eq. '41210' .or. 
     >    name(4:8) .eq. '41410' .or. 
     >    name(4:8) .eq. '41610' .or. 
     >    name(4:8) .eq. '42010' .or. 
     >    name(4:8) .eq. '42210' .or. 
     >    name(4:8) .eq. '42410' .or. 
     >    name(4:8) .eq. '42610' .or. 
     >    name(4:8) .eq. '42810' .or. 
     >    name(4:8) .eq. '43010' .or. 
     >    name(4:8) .eq. '43210' .or. 
     >    name(4:8) .eq. '43410' .or. 
     >    name(4:8) .eq. '50010' .or. 
     >    name(4:8) .eq. '50210' .or. 
     >    name(4:8) .eq. '50410' .or. 
     >    name(4:8) .eq. '50610' .or. 
     >    name(4:8) .eq. '50810' .or. 
     >    name(4:8) .eq. '51010' .or. 
     >    name(4:8) .eq. '51210' .or. 
     >    name(4:8) .eq. '51410' .or. 
     >    name(4:8) .eq. '51610' .or. 
     >    name(4:8) .eq. '51810' .or. 
     >    name(4:8) .eq. '52010' .or. 
     >    name(4:8) .eq. '52210' .or. 
     >    name(4:8) .eq. '52410' .or. 
     >    name(4:8) .eq. '52610' .or. 
     >    name(4:8) .eq. '52810' .or. 
     >    name(4:8) .eq. '53010' .or. 
     >    name(4:8) .eq. '53210' .or. 
     >    name(4:8) .eq. '53410' .or. 
     >    name(4:8) .eq. '60010' .or. 
     >    name(4:8) .eq. '60210' .or. 
     >    name(4:8) .eq. '60410' .or. 
     >    name(4:8) .eq. '60610' .or. 
     >    name(4:8) .eq. '60810' .or. 
     >    name(4:8) .eq. '61010' .or. 
     >    name(4:8) .eq. '61210' .or. 
     >    name(4:8) .eq. '61410' .or. 
     >    name(4:8) .eq. '61610' ) then
            name = 'QF1.'//name(4:8)
        endif   
      endif
      return
      end
