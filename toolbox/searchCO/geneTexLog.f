      character*120 cmmnd

C Create ./Log if does not exist already
      cmmnd = 'mkdir -p Log_CO'
      call system(cmmnd)

      cmmnd = 'cp zgoubi_searchCO-*.dat Log_CO'
      call system(cmmnd)

      cmmnd = 'cp  *eps ./Log_CO'
      write(*,*) '---------------------------------------'
      write(*,*) cmmnd
      call system(cmmnd)
 
      cmmnd = '~/zgoubi/struct/tools/searchCO/geneTexIni'
      write(*,*) '------------------------------------------'
      write(*,*) cmmnd
      call system(cmmnd)
      cmmnd = '~/zgoubi/struct/tools/searchCO/geneTexBody'
      write(*,*) '-------------------------------------------'
      write(*,*) cmmnd
      call system(cmmnd)
      cmmnd = '~/zgoubi/struct/tools/searchCO/geneTexEnd'
      write(*,*) '------------------------------------------'
      write(*,*) cmmnd
      call system(cmmnd)

      cmmnd = 
     >'cd Log_CO ; latex log_CO ; '
     >//'latex log_CO ; dvipdf log_CO'
      write(*,*) '-----------------------------------------------------'
      write(*,*) cmmnd
      call system(cmmnd)

      stop
      end
