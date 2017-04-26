      implicit double precision (a-h,o-z)
     
      parameter (lunR=11,lunW=12,lunDat=15,lunDa2=17,luntmp=14)
      character txt132*132, txt132c*132, let*1, txt4*4, txt4a*4, let1*1
      parameter (nCOmx=10001)
      dimension x(nCOmx),xp(nCOmx),z(nCOmx),
     >                   zp(nCOmx),s(nCOmx),d(nCOmx),let(nCOmx)
      character txtksy(10)*150, txt150*150
      character cmmnd*110
      character drctry*15, zgDatFile*50

      logical strcon
      INTEGER DEBSTR,FINSTR
      character typ2*12
      data typ2 / 'intrinsic' /
      parameter(pi=4.d0*atan(1.d0),c=2.99792458d8, deg2rd=pi/180.d0)
      parameter(zero=0.d0)
      data am, q, G / 938.27203d6,1.602176487d-19,1.7928474d0 /

      character(5) txtQz

C Input data ---------------------------------------------------------
      call system(
     >'cp -f geneZGDat4Xing.data geneZGDat4Xing.data_copy')
      open(unit=lunDat,file='geneZGDat4Xing.data')
      call readat(lunDat,
     >                   typ2,rJn2,circ,ah,gtr,Vp,phisD,am,q,G,ierr)
      write(*,fmt='(3a,9(1x,e12.4),1x,I4)') 
     >  ' Data read from geneZGDat4Xing.data : ',
     >  ' typ2, Jn2, circ, ah, gtr, Vp, phisD, am, q, G, ierr = ',
     >  typ2,rJn2,circ,ah,gtr,Vp,phisD,am,q,G,ierr
C typ2 = 'imperfection' or 'intrinsic' 
      close(lunDat)

C  geneZGDat4Xing.In is a copy of calcStrength.out 
      open(unit=lunDa2,file='geneZGDat4Xing.In')
C---------------------------------------------------------------------

      if(ierr.eq.1) then
        write(*,*) '  Failed to read a data from geneZGDat4Xing.data'
      endif
      if(ierr.ge.1) then 
        write(*,*) '  geneZGDat4Xing.data should exist and contain,', 
     >                                            ' in that order :'
        write(*,*) ' = M, algebraic_Qz :  such that resonant',
     >                                   ' gG=M + sign(Qz)|Qz|'
        write(*,*) ' ring circumference, harmoniq # '
        write(*,*) ' particle mass (MeV), charge (C), G'
        write(*,*) ' numb. of tracking turns upstream and downstream',
     >                                                  ' of resonance'
        stop
      endif
C-----------------------------
     
 20   continue
        read(lunDa2,fmt='(a)') txt132
        read(txt132,*) M
c        write(*,*) ' Press return to continue  M = ',M
c            read(*,*)
        read(txt132(5:12),*) txtQz
        read(txtQz,*,err=27,end=27) Qz
 27     continue
        Qz = 0.d0
 28     continue
c        write(*,*) ' Press return to continue   Qz= ',Qz 
c            read(*,*)
        read(txt132(12:132),*) aJn2im,aJ2in,ymu,zma,vkick
        write(*,*) ' Read from geneZGDat4Xing.In, ',
     >  'M,Qz,aJn2im,aJ2in,zma,vkick : ',M,Qz,aJn2im,aJ2in,zma,vkick
c        write(*,*) ' Press return to continue'
c            read(*,*)

        gamma = (dble(M) + Qz) / G
        if    (typ2 .eq. 'imperfection') then
C aJn2im is Jn^2/zmax^2
          aJn2 = aJn2im
          zma2 = 0.2d0 * rJn2 / aJn2
          if(zma2 .gt. 20.d-3) zma2=20.d-3
          vk = sqrt(zma2) * vkick/zma 
          epspi = zero
          dE = 300.d0 * sqrt(aJn2*zma2) / G * am
        elseif(typ2 .eq. 'intrinsic') then
C aJ2in is Nn^2/epsilon/pi
          aJn2 = aJ2in
          epspi = rJn2 / aJn2
          if(gamma*epspi .gt. 20.d-6) epspi = 20.d-6/gamma
          zma2 = 0.d0
          dE = 300.d0 * sqrt(aJn2*epspi) / G * am
        else
          write(*,*) 'typ2 = ',typ2,' :  No such choice available ! '
          stop 
        endif
C Note on value of nrbl : 
C Static,     Sz^2 = Delta^2/Jn^2 / (1 + Delta^2/Jn^2) 
C             Delta^2/Jn^2 = Sz^2 / (1 - Sz^2)
C  Distance to resoannce = gG-(n +/- Qz) = 7Jn corresponds to Sz=0.9900, i.e. 1% depolarization
C Sz=.9999 requires gG-(n +/- Qz) = 70.7Jn
C gG-(n +/- Qz) = 100Jn corresponds to Sz=0.9999500
C gG-(n +/- Qz) = 300Jn corresponds to Sz=0.99999444
C gG-(n +/- Qz) = 1000Jn corresponds to Sz=0.999999500

        E = gamma * am
        T = E - am
        p = sqrt(T * (T + 2.d0*am))
        boro = p/c
        if(gamma.gt.gtr) phisD = 180.d0 - phisD
        phis = phisD * deg2rd        
        dEturn = Vp * sin(phis)
        nrbl2 = nint(dE/dEturn)
        Ei = E - dble(nrbl2) * dEturn
        boroi = sqrt((Ei-am) * ((Ei-am) +2.d0*am))/c

        write(*,*) 
        write(*,*) ' M, Qz : ', M, Qz
        write(*,*) 
        write(*,*) ' Working conditions : '
        write(*,*) ' Mass, charge, G : ',am,q,G
        write(*,*) ' Emittance, max z_c.o. : ',epspi, sqrt(zma2)
        write(*,*) ' E, kin-E, Brho on resonance : ',E, T, boro
        write(*,*) ' Distance to the resonance : ',dE,' (eV)'
        write(*,*) ' E, Brho at start : ',Ei, boroi,', given  ',
     >    nrbl2,'  turns to resonance ' 
        write(*,*) ' Xing speed, GdE/2piM = ',G*dEturn/2.d0/pi/am
        write(*,*) 

        open(unit=lunR,file='zgoubi_geneZGDat4Xing-In.dat')

 2    continue

        rewind(lunR)
        open(unit=lunW,file='zgoubi_geneZGDat4Xing-Out.dat')
        rewind(lunW)

C Read in zgoubi_geneZGDat4Xing-In.dat
        read(lunR,fmt='(a)') txt132       ! title
        write(lunW,fmt='(a)') 
     >  'Data generated by  geneZGDat4Xing '
          read(lunR,fmt='(a)') txt132     ! OBJET
          write(lunW,*) txt132(debstr(txt132):finstr(txt132))
          read(lunR,fmt='(a)') txt132     ! BORO
          write(lunW,fmt='(f14.8,a2)') boroi,'d3'
          read(lunR,*) kobj
          write(*,*) ' KOBJ = ',kobj
          if(kobj .ne. 8) stop '** OBJET should use KOBJ=8 '
          write(lunW,*) '   8 ' 
          read(lunR,fmt='(a)') txt132
          write(lunW,*) txt132(debstr(txt132):finstr(txt132))
          read(lunR,fmt='(a)') txt132
          write(lunW,*) txt132(debstr(txt132):finstr(txt132))
          read(lunR,fmt='(a)') txt132
          write(lunW,*) txt132(debstr(txt132):finstr(txt132))
          read(lunR,*) alpz, btaz, dum
          write(lunW,fmt='(1p,3(1x,g14.6))') alpz, btaz, epspi
          read(lunR,fmt='(a)') txt132
          write(lunW,*) txt132(debstr(txt132):finstr(txt132))

      txt132 = '''PARTICUL'''
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))  
      write(lunW,fmt='(1p,3e17.9,a)') am/1d6, q , G, '  0. 0. '
      txt132 = '''SPNTRK'''
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))  
      txt132 = '  3       '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))  
      txt132 = '''FAISCEAU'''
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))  
      txt132 = '''PICKUPS'''
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))  
      txt132 = '   1        '
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))  
      txt132 = '   #Start'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))  
      txt132 = '''FAISTORE'''
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))  
      txt132 = '   b_zgoubi.fai  #End'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))  
      txt132 = '     1'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))  
      txt132 = '''MARKER''  #Start'
      write(lunW,*) txt132(debstr(txt132):finstr(txt132))  

C Complete zgoubi_geneZGDat4Xing-Out.dat from content of zgoubi_geneZGDat4Xing-In.dat
 1      continue
          read(lunR,fmt='(a)',end=10) txt132
          txt132 = txt132(debstr(txt132):132)

          if(strcon(txt132,'''PARTICUL''',10,
     >                                       IS)) then 
            read(lunR,fmt='(a)') txt132

          elseif(strcon(txt132,'''FAISCEAU''',10,
     >                                           IS)) then 

          elseif(strcon(txt132,'''SPNTRK''',8,
     >                                       IS)) then 
            read(lunR,fmt='(a)') txt132

          elseif(strcon(txt132,'''PICKUPS''',9,
     >                                       IS)) then 
            read(lunR,fmt='(a)') txt132
            read(lunR,fmt='(a)') txt132

          elseif(strcon(txt132,'''SPNSTORE''',10,
     >                                       IS)) then 
            read(lunR,fmt='(a)') txt132
            read(lunR,fmt='(a)') txt132

          elseif(strcon(txt132,'''FAISTORE''',10,
     >                                       IS)) then 
            read(lunR,fmt='(a)') txt132
            read(lunR,fmt='(a)') txt132

          elseif(strcon(txt132,'''SCALING''',9,
     >                                        IS)) then 
            write(lunW,*) txt132(debstr(txt132):finstr(txt132))
            read(lunR,fmt='(a)') txt132
            write(lunW,*) txt132(debstr(txt132):finstr(txt132))
            read(txt132,*) ii1, ii2
            do ii = 1, ii2
              read(lunR,fmt='(a)',end=10) txt132
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
              read(lunR,fmt='(a)',end=10) txt132
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
              read(lunR,fmt='(a)',end=10) txt132
              write(lunW,fmt='(f14.8)') boroi
c              write(*,fmt='(i4,f14.8)') ii2, boroi
c                   read(*,*) 
              read(lunR,fmt='(a)',end=10) txt132
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
            enddo

          elseif(strcon(txt132,'''MARKER''',8,
     >                                        IS)) then 
            txt132c = txt132(IS+8:132-8)
            txt132c = txt132c(debstr(txt132c):debstr(txt132c)+6)
            if(txt132c.eq.'#Start') then
            elseif(txt132c.eq.'#End') then
            else              
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
            endif

          elseif(strcon(txt132,'''TWISS''',7,
     >                                        IS)) then 
            read(lunR,fmt='(a)') txt132

          elseif(strcon(txt132,'DVCA02',8,                 ! AGS
     >                                        IS)) then 
            write(lunW,*) txt132(debstr(txt132):finstr(txt132))
            read(lunR,fmt='(a)',end=10) txt132
            write(lunW,*) txt132(debstr(txt132):finstr(txt132))
            read(lunR,*,end=10) xl,ro,b1
            b1 =  boroi * vk /xl        *1.d2 !*********
C             write(*,*) ' geneZGDat vk : ',vk,vkick,rJn2,aJn2,b1
C                   stop
            write(lunW,fmt='(1p,3(1x,g14.6),a)')
     >               xl,ro,b1,' 0. 0. 0. 0. 0. 0. 0. 0. 0.'

          elseif(strcon(txt132,'''CAVITE''',8,
     >                                        IS)) then
            read(lunR,fmt='(a)',end=10) txt132
            read(lunR,fmt='(a)',end=10) txt132
            read(lunR,fmt='(a)',end=10) txt132

          elseif(strcon(txt132,'''REBELOTE''',10,
     >                                         IS)) then
            read(lunR,fmt='(a)',end=10) txt132

          elseif(strcon(txt132,'''END''',5,
     >                                       IS)) then

            goto 11
          else
            
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
          endif

        goto 1

 11     continue
              txt132 = '''CAVITE'''
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
              txt132 = ' 2   .1'
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
              write(lunW,fmt='(f18.11,3x,f7.2)') circ, ah
              write(lunW,fmt='(1p,e16.8,3x,e20.12)') Vp, phis
              txt132 = '''MARKER''   #End'
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
              txt132 = '''REBELOTE'''
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
              write(lunW,fmt='(i8,a)') 2*nrbl2, '  0.2  99'

C When this zgoubi run finishes, Keyword 'SYSTEM' will cause launching of next zgoubi run
C Read file number from temporary storage by scanSpinResonances.f   
              open(unit=luntmp,file='scanSpinResonances.tmp')
              read(luntmp,*) ifile, drctry, zgDatFile
              close(luntmp)        
              txt132 = '''SYSTEM'''
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
              write(lunW,*) '  2'
              write(txt4,fmt='(i4)') ifile
              write(txt4a,fmt='(i4)') ifile+1
              cmmnd = 
     >        ' cd .. ; ' // 
     >        'sed -i ''s/'//txt4//' = xing simulation number/'
     >        //txt4a//' = xing simulation number/g'' '// 
     >        ' scanSpinResonances.count'
              write(*,*) cmmnd
              write(*,*) cmmnd
              write(*,*) cmmnd
              write(*,*) cmmnd
              write(lunW,*) cmmnd
              cmmnd = 'cd .. ; '//
     >        '~/zgoubi/toolbox/spin/resonanceXing/xing_scanResonances'
     >        //'/scanSpinResonances'
              write(lunW,*) cmmnd

              txt132 = '''END'''
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))

 10   continue
      write(*,*) ' '
      write(*,*) ' End of zgoubi_geneZGDat4Xing-In.dat input file has'
     >,' been reached, file zgoubi_geneZGDat4Xing-Out.dat completed.'
      write(*,*) ' ------------'
      write(*,*) ' '
      close(lunR)
      close(lunW)

C Build zgoubi_searchCO-In.dat and run searchCO
      open(unit=lunR,file='zgoubi_geneZGDat4Xing-Out.dat')
      open(unit=lunW,file='zgoubi_searchCO-In.dat')
C Replaces objet 8 by objet 2
          read(lunR,fmt='(a)',end=13) txt132   ! titl
          write(lunW,*) txt132(debstr(txt132):finstr(txt132))
          read(lunR,fmt='(a)',end=13) txt132   ! objet
          write(lunW,*) txt132(debstr(txt132):finstr(txt132))
          read(lunR,fmt='(a)',end=13) txt132   ! boro
          write(lunW,*) txt132(debstr(txt132):finstr(txt132))
          read(lunR,fmt='(a)',end=13) txt132   ! kobj
          write(lunW,*) ' 2 '
          read(lunR,fmt='(a)',end=13) txt132   !  1 1 1
          write(lunW,*) ' 1 1 '
          read(lunR,*,end=13) yo, ypo, zo, zpo, xo, do, let1
          write(lunW,fmt='(1p,4(e12.4,1x),e9.1,1x,e14.6,1x,3a1)') 
     >      yo*1d2, ypo*1d3, zo*1d2, zpo*1d3, xo, do, '''',let1,''''
          read(lunR,*,end=13)  alfy, bety, epsy
          write(lunW,*) ' 1 '
          read(lunR,*,end=13)  alfz, betz, epsz
          read(lunR,*,end=13)  alfl, betl, epsl
C Complete zgoubi_searchCO-In.dat from content of zgoubi_geneZGDat4Xing-In.dat
 12   continue
          read(lunR,fmt='(a)',end=13) txt132
          if(strcon(txt132,'''CAVITE''',8,
     >                                        IS)) then
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
              read(lunR,*,end=13) icav
              write(lunW,*) 0 
              read(lunR,fmt='(a)',end=13) txt132
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
              read(lunR,fmt='(a)',end=13) txt132
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
          elseif(strcon(txt132,'''REBELOTE''',10,
     >                                         IS)) then
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
              read(lunR,fmt='(a)',end=13) txt132
              write(lunW,*) ' 9  0.2  99'
          elseif(strcon(txt132,'''SYSTEM''',8,
     >                                        IS)) then
              read(lunR,*,end=13) ksys
              do kk = 1, ksys
                read(lunR,fmt='(a)',end=13) txtksy(kk)
              enddo
          else
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
          endif
          
        goto 12

 13     continue
      write(*,*) ' '
      write(*,*) ' End of zgoubi_geneZGDat4Xing-Out.dat input file has'
     >,' been reached, file zgoubi_searchCO-In.dat completed.'
      write(*,*) ' ------------'
      write(*,*) ' '
      close(lunR)
      close(lunW)

C Get closed orbit
      cmmnd = '~/zgoubi/toolbox/searchCO/searchCO_HV'
      call system(cmmnd)
      open(unit=lunR,file='zgoubi_searchCO-Out.dat')
      open(unit=lunW,file='zgoubi_geneZGDat4Xing-Out.dat')
C get co coordinates and rebuild objet in zgoubi_geneZGDat4Xing-Out.dat'
          read(lunR,fmt='(a)',end=16) txt132   ! titl
          write(lunW,*) txt132(debstr(txt132):finstr(txt132))
          read(lunR,fmt='(a)',end=16) txt132   ! objet
          write(lunW,*) txt132(debstr(txt132):finstr(txt132))
          read(lunR,fmt='(a)',end=16) txt132   ! boro
          write(lunW,*) txt132(debstr(txt132):finstr(txt132))
          read(lunR,fmt='(a)',end=16) txt132   ! kobj
          write(lunW,*) ' 8 '
          read(lunR,fmt='(a)',end=16) txt132   !  1 1 1
          write(lunW,*) ' 1 1 1 '
          read(lunR,*,end=16) yo, ypo, zo, zpo, xo, do, let1
          write(lunW,*) 
     >      yo/1d2, ypo/1d3, zo/1d2, zpo/1d3, xo, do, '''',let1,''''
          read(lunR,fmt='(a)',end=16) txt132   !  '  1  '
          write(lunW,*)  alfy, bety, epsy
          write(lunW,*)  alfz, betz, epsz
          write(lunW,*)  alfl, betl, epsl

 15       continue
            read(lunR,fmt='(a)',end=16) txt132
            if(strcon(txt132,'''CAVITE''',8,
     >                                        IS)) then
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
              read(lunR,fmt='(a)',end=16) txt132
              write(lunW,*) icav
              read(lunR,fmt='(a)',end=16) txt132
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
              read(lunR,fmt='(a)',end=16) txt132
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))

            elseif(strcon(txt132,'''REBELOTE''',10,
     >                                         IS)) then
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
              read(lunR,fmt='(a)',end=16) txt132
              write(lunW,*) 2*nrbl2, '  0.2  99'

              txt132 = '''SYSTEM'''
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
              write(lunW,*) ksys
              do kk = 1, ksys
                txt150 = txtksy(kk)
                write(lunW,fmt='(a)') 
     >             txt150(debstr(txt150):finstr(txt150))
              enddo

            else
              write(lunW,*) txt132(debstr(txt132):finstr(txt132))
            endif
            goto 15

 16     continue
      write(*,*) ' '
      write(*,*) ' End of zgoubi_searchCO-Out.dat input file has'
     >,' been reached, file zgoubi_geneZGDat4Xing-Out.dat completed.'
      write(*,*) ' ------------'
      write(*,*) ' '
      close(lunR)
      close(lunW)


      stop
      end
      FUNCTION STRCON(STR,STRIN,NCHAR,
     >                                IS)
      implicit double precision (a-h,o-z)
      LOGICAL STRCON
      CHARACTER STR*(*), STRIN*(*)
C     ------------------------------------------------------------------------
C     .TRUE. if the string STR contains the string STRIN with NCHAR characters
C     at least once.
C     IS = position of first occurence of STRIN in STR
C     ------------------------------------------------------------------------

      INTEGER DEBSTR,FINSTR

      II = 0
      DO 1 I = DEBSTR(STR), FINSTR(STR)
        II = II+1
        IF( STR(I:I+NCHAR-1) .EQ. STRIN ) THEN
          IS = II
          STRCON = .TRUE.
          RETURN
        ENDIF
 1    CONTINUE
      STRCON = .FALSE.
      RETURN
      END
      FUNCTION DEBSTR(STRING)
      implicit double precision (a-h,o-z)
      INTEGER DEBSTR
      CHARACTER * (*) STRING

C     --------------------------------------
C     RENVOIE DANS DEBSTR LE RANG DU
C     1-ER CHARACTER NON BLANC DE STRING,
C     OU BIEN 0 SI STRING EST VIDE ou BLANC.
C     --------------------------------------

      DEBSTR=0
      LENGTH=LEN(STRING)
C      LENGTH=LEN(STRING)+1
1     CONTINUE
        DEBSTR=DEBSTR+1
C        IF(DEBSTR .EQ. LENGTH) RETURN
C        IF (STRING(DEBSTR:DEBSTR) .EQ. ' ') GOTO 1
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
      implicit double precision (a-h,o-z)
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
      subroutine readat(lunDat,
     >                   typ2,rJn2,circ,ah,gtr,Vp,phisD,am,q,G,ierr)
      implicit double precision (a-h,o-z)
      character typ2*(*)
      ierr = 0
      read(lunDat,fmt='(a)',err=99,end=98) typ2
      read(lunDat,*,err=99,end=98) rJn2
      read(lunDat,*,err=99,end=98) circ, ah
      read(lunDat,*,err=99,end=98) gtr
      read(lunDat,*,err=99,end=98) Vp, phisD
      read(lunDat,*,err=99,end=98) am, q, G
      return
 99   continue
      ierr = 1
      write(*,*) ' error during read in data file'
      return
 98   continue
      write(*,*) ' End of data file reached'
      return
      end

