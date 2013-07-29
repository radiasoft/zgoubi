      character txt80*80, txtfig*20, txtM*4, txtQ*20
      parameter (lunW=8)
      integer debstr, finstr

      open(unit=lunW,file='./Log_CO/log_CO.tex',status='old')

      call gotoEnd(lunW,
     >                  line)

      txt80=''
      write(lunW,*) txt80
      txt80='\end{document}'
      write(lunW,*) txt80
      txt80=' '
      write(lunW,*) txt80
      txt80=' '
      write(lunW,*) txt80

      close(lunW)
      stop
      end
      subroutine gotoEnd(lunW,
     >                        line)
      character*132 txt132
      line = 0
 1    continue
        read(lunW,*,err=99,end=10) txt132
        line = line + 1
        goto 1

 10   continue
      backspace(lunW)
      return

 99   continue
      write(*,*) ' '
      write(*,*) '************* '
      write(*,*) 'Pgm geneTexEnd.  Error reading in lunW ' 
      write(*,*) 'line = ',line
      write(*,*) '************* '
      write(*,*) ' '
      stop 
      end

