C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C----- PLOT SPECTRUM     
      implicit double precision (a-h,o-z)
      CHARACTER txt132*132, txt9*9
      parameter (mxico=19999)
      dimension xco(MXICO),xpco(MXICO),zco(MXICO),zpco(MXICO)
      dimension dpco(MXICO),tof(MXICO)
      dimension path(MXICO), pp0(MXICO)
      dimension betx(MXICO), alfx(MXICO), Dx(MXICO), Dpx(MXICO)
      dimension betz(MXICO), alfz(MXICO), Dz(MXICO), Dpz(MXICO)
      dimension phix(MXICO), phiz(MXICO), numlm(MXICO)
      LOGICAL IDLUNI
      LOGICAL STRCON, EXS
      CHARACTER LET*1, txt1*1
      logical okpart
      CHARACTER txt35(MXICO)*35
      CHARACTER label*8, lbl1*8

      dimension YY(MXICO),TT(MXICO),ZZ(MXICO),PP(MXICO)

      INTEGER DEBSTR,FINSTR

      logical part
      logical gtText, ok

      data txt35 / MXICO*' ' /
      data label / 'all'/
      data okpart / .false. /

C EMMA
      data am, qe /  0.51099892, 1.602176487e-19  /

C Just a temp storage for passing K, Xi :
      INQUIRE(FILE='betaFromMatrix.in',exist=EXS)
      if(exs) then
        open(unit=34,file='betaFromMatrix.in')
        read(34,*,err=10,end=10)  label
        write(*,*) ' Filter label ***',label,'***'
        close(34)
      else
        label='all'
      endif
 10   continue

C Open zgoubi.res
      IF (IDLUNI(lunR)) THEN
        OPEN(UNIT=lunR,FILE='zgoubi.res',ERR=799)
      ELSE
        GOTO 798
      ENDIF

C Search zgoubi.res for co coordinates
c        read(lunR,fmt='(a)',end=62) txt132     ! title
c        read(lunR,fmt='(a)',end=62) txt132     ! 'OBJET'
        ok = gtText(lunR,'''OBJET''',
     >                               txt132)   
        read(lunR,*,end=62) BORO 
        read(lunR,fmt=*,end=62) objk
        kobj = objk
        kobj2 = nint(100.*objk) - 100*kobj
        if(kobj2.eq.0) then
          nCO = 1
        else
          nCO = kobj2
        endif
        read(lunR,fmt='(a)',end=62) txt132
        do ico = 1, nco
          read(lunR,*,end=62)
     >    xco(ico),xpco(ico),zco(ico),zpco(ico),dum,dpco(ico)
        enddo

C Search zgoubi.res for MATRIX outputs
      sDx = 0.d0
      sDz = 0.d0
      iico = 1
      jjco = 1
 61   continue

        if(iico.gt.MXICO) stop ' Too many elements. Increase MXICO'

        read(lunR,fmt='(a132)',end=62) txt132
c        write(*,fmt='(a132)') txt132

        if(strcon(txt132,'''PARTICUL''',10,
     >                                            IS)) then 
          okpart = .true.
          read(lunR,*,end=62) am, q

        elseif(strcon(txt132,'************************************',36,
     >                                                        IS)) then 
          read(lunR,fmt='(a132)',end=62) txt132
          txt35(iico) = txt132(29:37)
          read(txt132(38:finstr(txt132)),*,end=63) lbl1
C          write(*,*) txt1,lbl1,' *************'

          if(label(debstr(label):finstr(label)) .ne. 'all') then
c           write(*,*) lbl1(debstr(lbl1):finstr(lbl1)),
c     >                    label(debstr(label):finstr(label))
            if(lbl1(debstr(lbl1):finstr(lbl1)) .ne. 
     >         label(debstr(label):finstr(label))) then 

 64           read(lunR,fmt='(a132)',end=62) txt132
              if(.not.
     >        strcon(txt132,'************************************',36,
     >                                                     IS)) goto 64
              backspace(lunR)
              goto 61  
            endif
          endif

 63     continue

        elseif(strcon(txt132,'Reference particle (#',21,
     >                                                  IS)) then 
          txt132 = txt132(debstr(txt132):finstr(txt132))
          if(strcon(txt132,'), path length',14,
     >                                          IS)) then 
            read(txt132(44:56),*) path(iico) 
            read(txt132(84:96),*) pp0(iico) 
           
c              write(*,*) txt132(debstr(txt132):finstr(txt132))
c                write(*,*) iico, path(iico) ,pp0(iico)  
c                      read(*,*)
          endif

c        elseif(strcon(txt132,'Reference, absolute (part #',27,
        elseif(strcon(txt132,'Reference, before change of frame',33,
     >                                                    IS)) then 
            read(lunR,*) 
     >        D1,YY(iico),TT(iico),ZZ(iico),PP(iico),SS,time

c                write(*,*) 
c     >        D1,YY(iico),TT(iico),ZZ(iico),PP(iico),SS,time
c                     read(*,*)
           
        elseif(strcon(txt132,'BEAM  MATRIX (beta',18,
     >                                               IS)) then 

          if(strcon(txt132,'FINAL',5,
     >                                IS)) then 

c              write(*,*) '  Reading optical fuctions'
            read(lunR,fmt='(a)',end=62) txt132
            read(lunR,*,end=62) betx(iico),dum,dum,dum,dum,Dx(iico)
            read(lunR,*,end=62) alfx(iico),dum,dum,dum,dum,Dpx(iico)
c            alfx(iico) = -alfx(iico)
            read(lunR,*,end=62) dum,dum,betz(iico),dum,dum,Dz(iico)
            read(lunR,*,end=62) dum,dum,alfz(iico),dum,dum,Dpz(iico)
c            alfz(iico) = -alfz(iico)

            read(lunR,*,end=62) dum
            read(lunR,*,end=62) dum
            read(lunR,fmt='(a)',end=62) txt132
            read(lunR,fmt='(a)',end=62) txt132
            read(lunR,fmt='(a)',end=62) txt132

c              write(*,*) '  Reading betatron phase'
            read(lunR,*,end=62) phix(iico),phiz(iico)

            sDx = sDx + Dx(iico)
            sDz = sDz + Dz(iico)

            iico = iico + 1
          endif

        elseif(strcon(txt132,'TWISS  parameters,  periodicity  of',35,
     >                                                       IS)) then 

            read(lunR,fmt='(a)',end=62) txt132
            read(lunR,fmt='(a)',end=62) txt132
            read(lunR,fmt='(a)',end=62) txt132

            read(lunR,*,end=62) betx(iico),dum,dum,dum,dum,Dx(iico)
            read(lunR,*,end=62) alfx(iico),dum,dum,dum,dum,Dpx(iico)
            alfx(iico) = -alfx(iico)
            read(lunR,*,end=62) dum,dum,betz(iico),dum,dum,Dz(iico)
            read(lunR,*,end=62) dum,dum,alfz(iico),dum,dum,Dpz(iico)
            alfz(iico) = -alfz(iico)

            read(lunR,*,end=62) dum
            read(lunR,*,end=62) dum
            read(lunR,fmt='(a)',end=62) txt132
            read(lunR,fmt='(a)',end=62) txt132
            read(lunR,fmt='(a)',end=62) txt132

c              write(*,*) '  Reading betatron phase'
            read(lunR,fmt='(a)',end=62) txt132
c            write(*,*) txt132(27:37)
c            write(*,*) txt132(54:64)
            read(txt132(27:37),*,end=62,err=639) phix(iico)
            read(txt132(54:64),*,end=62,err=639) phiz(iico)

            sDx = sDx + Dx(iico)
            sDz = sDz + Dz(iico)

            iico = iico + 1
 639        continue

        endif

      goto 61
 62   close(lunR)

      call system('cat betaFromMatrix.out >> betaFromMatrix.out_old ; 
     >                   rm -f betaFromMatrix.out ')
            IF (IDLUNI(LUNW)) THEN
              OPEN(UNIT=LUNW,FILE='betaFromMatrix.out',ERR=699)
            ELSE
              GOTO 698
            ENDIF
      WRITE(LUNW,*) 
     >'% s, alfx, betx, Dx, Dpx, alfz, betz, Dz, Dpz, Sphx, Sphz, phx,
     >phz, p/p0, lmnt#'
      sDx2 = 0.d0
      sDz2 = 0.d0
      iico = iico-1 
      dsphx = 0.d0
      dsphz = 0.d0
      betxM = -1.d10
      betzM = -1.d10
      YYMa = -1.d10
      ZZMa = -1.d10
      DxMa = -1.d10
      DzMa = -1.d10
      YYMi = 1.d10
      ZZMi = 1.d10
      DxMi = 1.d10
      DzMi = 1.d10
      do ico = 1, iico
        ppp = BORO*0.299792458*dpco(ico)
        EKin = sqrt(ppp*ppp + AM*AM) - am
        if(ico.gt.1) then
          if(phix(ico) .lt. phix(ico-1)) dsphx = dsphx+1.d0
          if(phiz(ico) .lt. phiz(ico-1)) dsphz = dsphz+1.d0
        endif
        sphx = dsphx +  phix(ico)
        sphz = dsphz +  phiz(ico)
        if(debstr(txt35(ico)).le.0) txt35(ico)='*'
        WRITE(LUNW,179) path(ico), 
c     >  -alfx(ico),betx(ico),Dx(ico),Dpx(ico),
c     >  -alfz(ico),betz(ico),Dz(ico),Dpz(ico), 
     >  alfx(ico),betx(ico),Dx(ico),Dpx(ico),
     >  alfz(ico),betz(ico),Dz(ico),Dpz(ico), 
     >  sphx, sphz, phix(ico),phiz(ico), pp0(ico)
     >  ,txt35(ico)(debstr(txt35(ico)):finstr(txt35(ico)))
     >  ,YY(ico),TT(ico),ZZ(ico),PP(ico),ico
 179    FORMAT(1P,E14.6,13(1X,E13.5),1x,a,4(1X,E13.5),1X,i6)
        sDx2 = sDx2 + Dx(ico) * Dx(ico)
        sDz2 = sDz2 + Dz(ico) * Dz(ico)
        if(betx(ico).gt.betxM) betxM = betx(ico)
        if(betz(ico).gt.betzM) betzM = betz(ico)
        if(YY(ico).gt.YYMa) YYMa = YY(ico)
        if(ZZ(ico).gt.ZZMa) ZZMa = ZZ(ico)
        if(Dx(ico).gt.DxMa) DxMa = Dx(ico)
        if(Dz(ico).gt.DzMa) DzMa = Dz(ico)
        if(YY(ico).lt.YYMi) YYMi = YY(ico)
        if(ZZ(ico).lt.ZZMi) ZZMi = ZZ(ico)
        if(Dx(ico).lt.DxMi) DxMi = Dx(ico)
        if(Dz(ico).lt.DzMi) DzMi = Dz(ico)
      enddo

      sDx = sDx / dble(iico)
      sDz = sDz / dble(iico)
      sDx2 = sDx2 / dble(iico)
      sDz2 = sDz2 / dble(iico)

      sigDx =  sqrt(sDx2 - sDx*sDx)
      sigDz =  sqrt(sDz2 - sDz*sDz)

      write(lunW,fmt='(2a)') 
     > ' averages : Dx, Dz, sigDx, sigDz, sqrt(sDx2), sqrt(sDz2) '
      write(lunW,fmt='(a,1p,6(1x,e14.6),1(1x,i4))') 
     > ' # ', sDx, sDz, sigDx, sigDz, sqrt(sDx2), sqrt(sDz2), iico
      write(lunW,fmt='(2a)') 
     > ' extrema : max. betx, betz, Dx, Dz ; min. Dx, Dz'
     > ,' ;  max.  xco, zco ; min.  xco, zco'
      write(lunW,fmt='(a,1p,10(1x,e14.6))') 
     > ' # ', betxM, betzM, DxMa, DzMa, DxMi, DzMi
     > ,YYMa, ZZMa, YYMi, ZZMi


      write(6,*) ' Job ended.  Went on well it seems...'
      goto 99

 698  WRITE(6,*) ' *** Problem : No idle unit for betaFromMatrix.out '
      GOTO 99
 699  WRITE(6,*) ' *** Problem at OPEN betaFromMatrix.out '
      GOTO 99

 798  WRITE(6,*) ' *** Problem : No idle unit for open zgoubi.res '
      GOTO 99
 799  WRITE(6,*) ' *** Problem at OPEN zgoubi.res '
      GOTO 99

 99   continue

      stop
      end
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

      function gtText(lunR,txt,
     >                         txt132)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical gtText
      character(*) txt
      character(*) txt132
      logical strcon
      integer debstr, finstr

      nt = len(txt(debstr(txt):finstr(txt)))

      read(lunR,fmt='(a)',err=99,end=98)  txt132

      dowhile(.not.strcon(txt132,txt(debstr(txt):finstr(txt)),nt, 
     >                                  IS))
        read(lunR,fmt='(a)',err=99,end=98)  txt132
      enddo

      gtText = .true.
      goto 10

 99   continue
      write(*,*) ' PGM searchCO_only_0000/gtText - 99 upon read.' 
      write(*,*) ' text was : ',txt(debstr(txt):finstr(txt))
      gtText = .false.
      goto 10

 98   continue
      write(*,*) ' PGM searchCO_only_0000/gtText - EOF upon read.' 
      write(*,*) ' text was : ',txt(debstr(txt):finstr(txt))
      gtText = .false.
      goto 10

 10   continue
      return
      end
