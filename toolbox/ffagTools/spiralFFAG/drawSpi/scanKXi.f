      logical exs, strcon
      character txt132*132, txt10*10
      character rep*1, prtcl*6
      data lunIn, lunOut, lunRes, lunData / 7, 8, 11, 13 / 
      data prtcl / 'proton' / 

      write(6,*) 
      write(6,*) '--------------------------------------------'
      write(6,*) ' scanKXi procedure now going on...' 

C---------- Input data ------------------
      open(unit=lunData,file='scanKXi.data')
      read(lunData,*) xKmi,dK,nK
      close(lunData)
C----------------------------------------

C Create a (new) output file 
      inquire(file='scanKXi.out',EXIST=EXS)
c      pause
      if (exs) then
        open(unit=lunOut,name='scanKXi.out')
        close(unit=lunOut,status='delete')
      endif
c      pause
      open(unit=lunOut,name='scanKXi.out')
c      pause
      write(lunOut,*) '   '

c      dK = (xKma - xKmi)/float(nK-1)
      xKi = xKmi - dK
      do k = 1, nK
        xK = xKi + float(k) *dK

        write(6,*)
        write(6,*) ' K = ',k
        write(6,*) ' '
 
C Change K value in lattice.data
        call system('cp -f lattice.data lattice.data_temp')
        call system('cp -f lattice.data dum')
        call change('lattice.data','lattice.data_temp',xK,
     >                             xiDegI,nCellI,r0I,pfI,gapI,prtcl)

C Run MAD  
        call system('rm -f dict ; ln -s 
     >   ~/mad/source/mad8.Linux/mad8.dict ./dict ;
     >   ~/mad/source/mad8.Linux/mad8 < lattice.cmd')

C Open input file 
        open(unit=lunIn,name='print')

        write(lunOut,*) ' ' 
        write(lunOut,fmt='(''nCell'',T10,''K'',T20,''xi'',T30,''Qx''
     >  ,T40,''Qy'',T47,''MaxbetX'',T57,''MaxBety'',T68,''MaxDx'')') 

        iter = 1
 10     continue
C Looking for 'Element            attribute'. That's where parameters can be read
        read(lunIn,fmt='(A)',end=98,err=99) txt132

        if    (strcon(txt132,'Element            attribute',28,
     >                                                       IS)) then

          read(lunIn,fmt='(A)',end=98,err=99) txt132
          read(txt132,fmt='(A40,E16.4)',end=98,err=99) txt10, cellN
          read(lunIn,fmt='(A)',end=98,err=99) txt132
          read(txt132,fmt='(A40,E16.4)',end=98,err=99) txt10, aK
          read(lunIn,fmt='(A)',end=98,err=99) txt132
          read(txt132,fmt='(A40,E16.4)',end=98,err=99) txt10, xiDeg

          do ir = 1, 11
            read(lunIn,fmt='(A)',end=98,err=99) txt132
          enddo

          read(lunIn,fmt='(A)',end=98,err=99) txt132

          if    (strcon(txt132,'total length =',14,
     >                                             IS)) then

            read(txt132,fmt='(A21,E12.4,A30,E12.4,A30,E12.4)',
     >         end=98,err=99) txt10,xLen,txt10,Qx,txt10,Qy
            read(lunIn,fmt='(A)',end=98,err=99) txt132
            read(lunIn,fmt='(A)',end=98,err=99) txt132
            read(txt132,fmt='(A21,E12.4,A30,E12.4,A30,E12.4)',
     >         end=98,err=99) txt10,alpha,txt10,betxM,txt10,betyM
            read(lunIn,fmt='(A)',end=98,err=99) txt132
            read(txt132,fmt='(A21,E12.4,A30,E12.4,A30,E12.4)',
     >         end=98,err=99) txt10,gtr,txt10,DxM,txt10,DyM
          
            nCell = nint(cellN)
            write(lunOut,fmt='(i3, 7f10.5, 2i4)') 
CC     >      nCell, aK, xiDeg, Qx/nCell, Qy/nCell, betxM, betyM, 
     >        nCell, aK, xiDeg, Qx, Qy, betxM, betyM,DxM,iter,K
             
          endif
       
          iter = iter + 1
        endif

        goto 10


 98     write(6,*) '  * End of job upon EOFile *    '
c        goto 999
        
 99     write(6,*) '  * End of job upon read error *   '
c        goto 999

        write(lunOut,*) '  ' 
        close(unit=lunIn)
      enddo

 999  continue

C Run MAD for nominal K, xi, will produce mad.ps for geneTex
      call system('rm -f dict ; ln -s 
     > ~/mad/source/mad8.Linux/mad8.dict ./dict ;
     > ~/mad/source/mad8.Linux/mad8 < lattice.cmd-givenKxi')
      call prNomi(lunIn, lunOut)

      close(unit=lunOut)


        call system('gnuplot 
     >~/mad/structure/ffag/tools/spiralFFAG/scanKXi/gnuplotBtaXMax.cmd')
c        pause
        call system('gnuplot 
     >~/mad/structure/ffag/tools/spiralFFAG/scanKXi/gnuplot3DBtaXMa.cmd'
     >)
c        pause
        call system('gnuplot 
     > ~/mad/structure/ffag/tools/spiralFFAG/scanKXi/gnuplotBtaZMax.cmd'
     >)
c        pause
        call system('gnuplot 
     >~/mad/structure/ffag/tools/spiralFFAG/scanKXi/gnuplot3DBtaZMa.cmd'
     >)
c        pause
        call system('gnuplot 
     > ~/mad/structure/ffag/tools/spiralFFAG/scanKXi/gnuplotBtaMiMa.cmd'
     >)
c        pause
        call system('gnuplot 
     > ~/mad/structure/ffag/tools/spiralFFAG/scanKXi/gnuplotQxQz.cmd')

        if(prtcl .eq. 'carbon') then 
          call system('gnuplot 
     >  ~/mad/structure/ffag/tools/spiralFFAG/scanKXi/gnuplotKxi_C.cmd')
        else
          call system('gnuplot 
     >    ~/mad/structure/ffag/tools/spiralFFAG/scanKXi/gnuplotKxi.cmd')
        endif        

c        pause
c        call system('gv mad.ps &')
c        call system('gv gnuplotBtaMax.eps &')
c        call system('gv gnuplotBtaMiMa.eps &')
c        call system('gv gnuplotQxQz.eps &')
c        call system('gv gnuplotKxi.eps &')

c        call system('rm lattice.cmd_temp')

      stop
      end

      FUNCTION STRCON(STR,STRIN,NCHAR,
     >                                IS)
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
      subroutine change(nameIn,nameOut,xK, 
     >                           xiDeg,nCell,r0,pf,gap,prtcl)
      character*(*) nameIn, nameOut, prtcl
      character*132 txt132
      logical strcon
      integer debstr
      open(unit=9,name=nameIn)
      open(unit=10,name=nameOut)
      read(9,fmt='(a)',end=99)  txt132  !! reads comment line. Next line is to be K
      write(10,fmt='(a)') txt132
      read(9,fmt='(a)',end=99) txt132   !! line that contains K=***
      write(10,fmt='(a,1p,g12.4,a)') 'K = ',xK, '!do not move this line'
      read(9,fmt='(a)',end=99)  txt132
      write(10,fmt='(a)') txt132
      if (STRCON(txt132,'=',1,
     >                       IS)) then
        read(txt132(IS+1:132),*) xiDeg
      endif
      read(9,fmt='(a)',end=99)  txt132
      write(10,fmt='(a)') txt132
      write(6,fmt='(a)') txt132
      if (STRCON(txt132,'=',1,
     >                       IS)) then
        read(txt132(IS+1:132),*) nCell
      endif 
      read(9,fmt='(a)',end=99)  txt132
      write(10,fmt='(a)') txt132
      write(6,fmt='(a)') txt132
      if (STRCON(txt132,'=',1,
     >                       IS)) then
        read(txt132(IS+1:132),*) r0
      endif
      read(9,fmt='(a)',end=99)  txt132
      write(10,fmt='(a)') txt132
      if (STRCON(txt132,'=',1,
     >                       IS)) then
        read(txt132(IS+1:132),*) pf
      endif
      read(9,fmt='(a)',end=99)  txt132
      write(10,fmt='(a)') txt132
      if (STRCON(txt132,'=',1,
     >                       IS)) then
        read(txt132(IS+1:132),*) gap
      endif
 1    read(9,fmt='(a)',end=99)  txt132
      write(10,fmt='(a)') txt132
      write(*,fmt='(a)') txt132
      if    (strcon(txt132,'carbon',6,
     >                                IS)) then
        prtcl = 'carbon'
      endif
      goto 1
 99   close(9)
      close(10)
        
      return
      end
      subroutine prNomi(lunIn, lunOut)
      logical exs, strcon
      character txt132*132, txt10*10
C Open input file 
      open(unit=lunIn,name='print')
      open(unit=lunOut,name='scanKXi.out')
      call go2end(lunOut)
        write(lunOut,*) ' ' 
        write(lunOut,fmt='(''nCell'',T10,''K'',T20,''xi'',T30,''Qx''
     >  ,T40,''Qy'',T47,''MaxbetX'',T57,''MaxBety'',T68,''MaxDx'')') 

        iter = 1
 10     continue
C Looking for 'Element            attribute'. That's where parameters can be read
        read(lunIn,fmt='(A)',end=98,err=99) txt132
        if    (strcon(txt132,'Element            attribute',28,
     >                                                       IS)) then

          read(lunIn,fmt='(A)',end=98,err=99) txt132
          read(txt132,fmt='(A40,E16.4)',end=98,err=99) txt10, cellN
          read(lunIn,fmt='(A)',end=98,err=99) txt132
          read(txt132,fmt='(A40,E16.4)',end=98,err=99) txt10, aK
          read(lunIn,fmt='(A)',end=98,err=99) txt132
          read(txt132,fmt='(A40,E16.4)',end=98,err=99) txt10, xiDeg

          do ir = 1, 11
            read(lunIn,fmt='(A)',end=98,err=99) txt132
          enddo

          read(lunIn,fmt='(A)',end=98,err=99) txt132

          if    (strcon(txt132,'total length =',14,
     >                                             IS)) then

            read(txt132,fmt='(A21,E12.4,A30,E12.4,A30,E12.4)',
     >         end=98,err=99) txt10,xLen,txt10,Qx,txt10,Qy
            read(lunIn,fmt='(A)',end=98,err=99) txt132
            read(lunIn,fmt='(A)',end=98,err=99) txt132
            read(txt132,fmt='(A21,E12.4,A30,E12.4,A30,E12.4)',
     >         end=98,err=99) txt10,alpha,txt10,betxM,txt10,betyM
            read(lunIn,fmt='(A)',end=98,err=99) txt132
            read(txt132,fmt='(A21,E12.4,A30,E12.4,A30,E12.4)',
     >         end=98,err=99) txt10,gtr,txt10,DxM,txt10,DyM
          
            nCell = nint(cellN)
            write(lunOut,fmt='(i3, 7f10.5,a)') 
CC     >      nCell, aK, xiDeg, Qx/nCell, Qy/nCell, betxM, betyM, 
     >      nCell, aK, xiDeg, nCell*Qx, nCell*Qy, betxM, betyM,DxM,
     >         ' 999 999'

            write(*,*) '  scanKXi, prNomi : '
            write(*,fmt='(i3, 7f10.5, a)') 
CC     >      nCell, aK, xiDeg, Qx, Qy, betxM, betyM,DxM,' 999 999'
     >      nCell, aK, xiDeg, nCell*Qx, nCell*Qy, betxM, betyM,DxM,
     >         ' 999 999'

          endif
       
          iter = iter + 1
        endif

        goto 10


 98     write(6,*) '  * End of job upon EOFile *    '
        goto 999
        
 99     write(6,*) '  * End of job upon read error *   '
        goto 999

c        write(lunOut,*) '  ' 
        close(unit=lunIn)

 999  continue
      return
      end
      subroutine go2end(lun)
      implicit double precision (a-h, o-z)
      character*132 txt
 1    continue
        read(lun,*,end=99,err=99) txt
        goto 1
 99   return
      end
