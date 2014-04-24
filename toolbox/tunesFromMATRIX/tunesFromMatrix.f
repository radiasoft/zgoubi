C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
      implicit double precision (a-h,o-z)
      CHARACTER txt132*132, txt9*9
      dimension xco(999), xpco(999), zco(999), zpco(999), 
     >                                  dpco(999), tof(999)
      dimension path(999)
      dimension betx(999), alfx(999), Dx(999), Dpx(999)
      dimension betz(999), alfz(999), Dz(999), Dpz(999)
      dimension xnu(999), znu(999)
      LOGICAL IDLUNI
      LOGICAL STRCON, EXS
      CHARACTER LET*1
      CHARACTER*14 txtxnu, txtznu
      logical okpart
      data okpart / .false. /

      logical part
C EMMA
      data am, qe /  0.51099892, 1.602176487e-19  /

C Just a temp storage for passing K, Xi :
      INQUIRE(FILE='tempKXi.dum',exist=EXS)
      if(exs) then
        open(unit=34,file='tempKXi.dum')
        read(34,*,err=10,end=10)  xK, xiDeg
        close(34)
      else
        xK= 0.
        xiDeg = 0. 
      endif
 10   continue

C Open zgoubi.res
      IF (IDLUNI(lunR)) THEN
        OPEN(UNIT=lunR,FILE='zgoubi.res',ERR=799)
      ELSE
        GOTO 798
      ENDIF

C Search zgoubi.res for orbit coordinates
        read(lunR,fmt='(a)',end=62) txt132     ! title
        read(lunR,fmt='(a)',end=62) txt132     ! 'OBJET'
        read(lunR,*,end=62) BORO 
        read(lunR,fmt='(a)',end=62) txt132     ! KOBJ.KOBJ2
        read(txt132,fmt=*,end=62) objk
        kobj = objk
        if(strcon(txt132,'.',1,
     >                         IS)) then
          read(txt132(is+1:132),*) kobj2
        else
          kobj2 = 0
        endif         
c        kobj = objk
c        kobj2 = nint(100.*objk) - 100*kobj
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
      iico = 1
      jjco = 1
      okpart = .false.
 61   continue

        read(lunR,fmt='(a132)',end=62) txt132

        if(.not. okpart) then
          if(strcon(txt132,'''PARTICUL''',10,
     >                                            IS)) then 
            okpart = .true.
            read(lunR,*,end=62) am, q
          endif
        endif

        if(strcon(txt132,'Beam  matrix  (beta',19,
     >                                            IS)) then 
          read(lunR,fmt='(a)',end=62) txt132
          read(lunR,*,end=62) betx(iico),dum,dum,dum,dum,Dx(iico)
          read(lunR,*,end=62) alfx(iico),dum,dum,dum,dum,Dpx(iico)
          alfx(iico) = -alfx(iico)
          read(lunR,*,end=62) dum,dum,betz(iico),dum,dum,Dz(iico)
          read(lunR,*,end=62) dum,dum,alfz(iico),dum,dum,Dpz(iico)
          alfz(iico) = -alfz(iico)
          read(lunR,fmt='(a)',end=62) txt132
          read(lunR,fmt='(a)',end=62) txt132
          read(lunR,fmt='(a)',end=62) txt132
          read(lunR,fmt='(a)',end=62) txt132
          read(lunR,fmt='(a)',end=62) txt132
          read(lunR,*,end=62) 
     >          txt9,txt9, txtxnu, txt9,txt9, txtznu
          if(strcon(txtxnu,'undefined',9,
     >                                   IS)) then 
            xnu(iico) = -1.d0
          else
            read(txtxnu,*)  xnu(iico)
          endif
          if(strcon(txtznu,'undefined',9,
     >                                   IS)) then 
            znu(iico) = -1.d0            
          else
            read(txtznu,*)  znu(iico)
          endif

C Test ----------------------------------------------------
            write(*,*) ' QY, QZ : ',xnu(iico),znu(iico),iico,nco
C End test ----------------------------------------------------

          iico = iico + 1
        endif

        if(iico .gt. 1) then      !! means MATRIX has already been seen
C This requires 'FAISCEAU' following 'MATRIX'
          if(strcon(txt132,'Time of flight',14,
     >                                         IS)) then 

c  read path length of traj 1 (modulo 11)
            backspace(lunR)
            backspace(lunR)
            read(lunR,fmt='(a)',end=62) txt132
            read(txt132,*) LET,IEX,fo1,fo2,fo3,fo4,fo5,fo6
     >      ,F1,F2,F3,F4,F5,F6
c            read(txt132,101) LET,IEX,fo1,fo2,fo3,fo4,fo5,fo6
c     >      ,F1,F2,F3,F4,F5,F6,  I
 101    FORMAT(A1,1X,I3,1X,F8.4,5F10.3,5X,F8.4,4F9.3,1X,F11.3,1X,I6)
            path(jjco) = f6
c  read tof of traj 1 (modulo 11)
            read(lunR,fmt='(a)',end=62) txt132
            read(txt132(39:132),*) tof(jjco)

c  jump to traj. 10 (modulo 11)
            do jj = 1, 16
              read(lunR,fmt='(a)',end=62) txt132
            enddo

c  read path length of traj 10 and 11 (modulo 11)
            read(lunR,fmt='(a)',end=62) txt132
C            read(txt132,101) LET,IEX,fo1,fo2,fo3,fo4,fo5,fo6
C     >      ,F1p,F2,F3,F4,F5,F6p,  I
            read(txt132,*) LET,IEX,fo1,fo2,fo3,fo4,fo5,fo6
     >      ,F1p,F2,F3,F4,F5,F6p,  I
c            write(*,*) ' Path length of traj. # ',I,' is ',F6p,' cm'
            read(lunR,fmt='(a)',end=62) txt132
            read(txt132(39:132),*) tofp
            read(lunR,fmt='(a)',end=62) txt132
C            read(txt132,101) LET,IEX,fo1,fo2,fo3,fo4,fo5,fo6
C     >      ,F1m,F2,F3,F4,F5,F6m,  I
            read(txt132,*) LET,IEX,fo1,fo2,fo3,fo4,fo5,fo6
     >      ,F1m,F2,F3,F4,F5,F6m,  I
c            write(*,*) ' Path length of traj. # ',I,' is ',F6m,' cm'
            read(lunR,fmt='(a)',end=62) txt132
            read(txt132(39:132),*) tofm

            jjco = jjco + 1
          endif
        endif

      goto 61
 62   close(lunR)

      call system('cat tunesFromMatrix.out >> tunesFromMatrix.out_old ; 
     >                   rm -f tunesFromMatrix.out ')
            IF (IDLUNI(LUNW)) THEN
              OPEN(UNIT=LUNW,FILE='tunesFromMatrix.out',ERR=699)
            ELSE
              GOTO 698
            ENDIF
            IF (IDLUNI(LUNW2)) THEN
              OPEN(UNIT=LUNW2,FILE='tunesFromMatrix.out2',ERR=699)
            ELSE
              GOTO 698
            ENDIF
      WRITE(LUNW,fmt='(a)') 
     >'#   1        2       3   4    5     6  '
     >//'  7     8     9     10     '
     >// '    11        12          13    '
      WRITE(LUNW,fmt='(a)') 
     >'# xco[m], x''co[rd], Qx, Qy, alfx, betx,'
     >//' alfy, bety, Dx, tof[mu_s], '
     >//'path[cm], E_Kin[MeV], p/p0(co)'
      WRITE(LUNW,fmt='(a)') '# '
      WRITE(LUNW2,fmt='(a)') 'E_kin(eV) xco(m) xpco(rad) TOF(s) Qx Qy'

      do ico = 1, nco
        pp = BORO*0.299792458*dpco(ico)
        EKin = sqrt(pp*pp + AM*AM) - am
        WRITE(LUNW,179) xco(ico)/1.d2, xpco(ico)/1.d3, 
     >  XNU(ico), ZNU(ico), alfx(ico), betx(ico), alfz(ico), betz(ico),
     >  Dx(ico), tof(ico), path(ico), EKin, dpco(ico) !!!!, xK, xiDeg
 179    FORMAT(1P,13(E16.7,1X))
        WRITE(LUNW2,180) EKin*1d6, xco(ico)*1.d-2, xpco(ico)*1.d-3, 
     >  tof(ico)*1d-6, XNU(ico),ZNU(ico)
 180    FORMAT(1P,6G16.8)
      enddo
      close(LUNW)
      close(LUNW2)

C Long write-up
      call system('cat tunesFromMatrix_LW.out >> 
     >tunesFromMatrix_LW.out_old ; 
     >rm -f tunesFromMatrix_LW.out ')
            IF (IDLUNI(LUNW)) THEN
              OPEN(UNIT=LUNW,FILE='tunesFromMatrix_LW.out',ERR=699)
            ELSE
              GOTO 698
            ENDIF
      WRITE(LUNW,fmt='(a)') 
     >'# xco, xpco, zco, zpco, path(s), dpco, tof(s), '//
     >' alfx, betx, alfz, betz, alfl, betl, '//
     >' Dx, Dpx, Dy, Dpy, Dl, Dpl, xmu, zmu, smu, E_Kin'

      alfl = 0.d0
      betl = 0.d0
      Dl = 0.d0
      Dpl = 0.d0
      smu = 0.d0
      do ico = 1, nco
        pp = BORO*0.299792458*dpco(ico)
        EKin = sqrt(pp*pp + AM*AM) - am
        WRITE(LUNW,149) xco(ico)/1.d2, xpco(ico)/1.d3, 
     >  zco(ico)/1.d2, zpco(ico)/1.d3, path(ico)/1d2,dpco(ico),tof(ico), 
     >  alfx(ico), betx(ico), alfz(ico), betz(ico),alfl, betl,
     >  Dx(ico), Dpx(ico), Dz(ico), Dpz(ico), Dl, Dpl, 
     >  XNU(ico), ZNU(ico), smu, EKin
 149    FORMAT(1P,23(E16.8,1X))
      enddo
      close(LUNW)


      write(6,*) ' Job ended.  Went on well it seems...'
      goto 99

 698  WRITE(6,*) ' *** Problem : No idle unit for tunesFromMatrix.out '
      GOTO 99
 699  WRITE(6,*) ' *** Problem at OPEN tunesFromMatrix.out '
      GOTO 99

 798  WRITE(6,*) ' *** Problem : No idle unit for open zgoubi.res '
      GOTO 99
 799  WRITE(6,*) ' *** Problem at OPEN zgoubi.res '
      GOTO 99

 99   continue
      if(.not. okpart) then
          write(*,*) ' '
          write(*,*) ' WARNING : no ''PARTICUL'' in zgoubi.dat'
          write(*,*) '           TOF won''t be computed '
          write(*,*) ' '
      endif


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
