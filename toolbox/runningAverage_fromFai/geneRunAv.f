      character(4) txRun

      open(2,file='runAv_1024.out')

      do i = 0, 1023
        write(txRun,fmt='(i0)') i 
        write(2,*) 'echo " Run'//trim(txRun)//'"'
        write(2,*) '( cd Run'//trim(txRun)
     >  //' ; ../runningAverage_fromFai )'
      enddo

         stop
             end

