      logical strcon, ok
      character(132) txt132
      character(*) fname
      character(*) titl
      integer debstr

      ierr = 0
      
      open(unit=lunIn,file=fname,status='old',err=98)

C Default mass number and charge number
      A = 1.d0
      Q = 1.d0

      read(lunIn,fmt='(a)',end=96,err=97)  txt132  !! reads comment line
      titl = txt132
      write(*,*) titl
 
C For compatibility with earlier version of dataFile
      T0 = -1.d0

 1    continue

      read(lunIn,fmt='(a)',end=96,err=97)  txt132
      ideb = debstr(txt132)
      txt132 = txt132(ideb:132)

      if (STRCON(txt132(1:6),'K',1,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) aK
        write(*,*) ' aK = ', aK
      elseif (STRCON(txt132(1:6),'xiDeg',5,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) ZTA
        write(*,*) ' xiDeg = ', ZTA
      elseif (STRCON(txt132(1:6),'nCell',5,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) nCell
        write(*,*) ' nCell = ', nCell
      elseif (STRCON(txt132(1:6),'r0',2,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) r0
        write(*,*) ' r0 = ', r0
      elseif (STRCON(txt132(1:6),'pf',2,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) pf
        write(*,*) ' pf = ', pf
      elseif (STRCON(txt132(1:6),'gap',3,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) gap
        write(*,*) ' gap = ', gap
      elseif (STRCON(txt132(1:6),'kappa',5,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) akappa
        write(*,*) ' kappa = ', akappa
      elseif (STRCON(txt132(1:6),'AA',2,
     >                                 IS)) then
        ok=STRCON(txt132,'=',1,
     >                         IS)
        read(txt132(IS+1:132),*) AA
        write(*,*) ' AA = ', A
      elseif (STRCON(txt132(1:6),'QQ',2,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) QQ
        write(*,*) ' QQ = ', A
      elseif (STRCON(txt132(1:6),'T1',2,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) T1
        write(*,*) ' T1 = ', T1
      elseif (STRCON(txt132(1:6),'T2',2,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) T2
        write(*,*) ' T2 = ', T2
      elseif (STRCON(txt132(1:6),'T0',2,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) T0
        write(*,*) ' T0 = ', T0
      elseif (STRCON(txt132(1:6),'nCO',3,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) nCO
        write(*,*) ' nCO = ', nCO
      elseif (STRCON(txt132(1:6),'map',3,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) map
        write(*,*) ' map = ', map
      elseif (STRCON(txt132(1:6),'kaseV',5,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) kaseV
        write(*,*) ' kaseV = ', kaseV
      elseif (STRCON(txt132(1:6),'precX',5,
     >                       IS)
     >  .and. (.not. STRCON(txt132(1:7),'precXp',6,
     >                                      IS))) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) precX
        write(*,*) ' precX = ', precX
		    read(*,*)
      elseif (STRCON(txt132(1:7),'precXp',6,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) precXp
        write(*,*) ' precXp = ', precXp
				    read(*,*)
      elseif (STRCON(txt132(1:6),'DltACN',6,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) DltACN
        write(*,*) ' DltACN = ', DltACN
      elseif (STRCON(txt132(1:6),'AT',2,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) AT
        write(*,*) ' AT = ', AT
      elseif (STRCON(txt132(1:6),'kDrift',6,
     >                       IS)) then
        ok=STRCON(txt132,'=',1,
     >                       IS)
        read(txt132(IS+1:132),*) kDrift
        write(*,*) ' kDrift = ', kDrift
      endif

      goto 1

 96   continue
      ierr = 1
      close(lunIn)
      write(*,*) ' readat ended upon read/eof'
      goto 99

 97   continue
      ierr = 2
      close(lunIn)
      write(*,*) ' readat ended upon read/error'
      goto 99

 98   continue
      ierr = 3
      close(lunIn)
      write(*,*) ' readat ended upon file open/error'
      goto 99

 99   continue
C For compatibility with earlier version of data file 
      if(T0 .le. 0.d0) T0 = T2


