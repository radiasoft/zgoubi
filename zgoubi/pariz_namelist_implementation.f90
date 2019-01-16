submodule(pariz_namelist_interface) pariz_namelist_implementation
  use assertions_interface, only : assert
  implicit none

contains

  module procedure initialize_input_parameters

    namelist/pariz/ ID, IZ, MMAP, MXX, MXY

    block 
      integer, parameter :: success=0
      integer file_unit, io_status
      logical pariz_nml_open

      inquire(file=pariz_namelist_file, number=file_unit, opened=pariz_nml_open)
      call assert( .not. pariz_nml_open, "read_from_file: file not already open" )
        !! If the above assertion is removed, then add code the following form:
        !! if(pariz_nml_open) then; rewind(...); else; open(...)
    
      open(newunit=file_unit,file=pariz_namelist_file,iostat=io_status,status='old')
      if (io_status/=success) return
        !! return with the default values if the namelist file is not found
  
      read(file_unit,nml=pariz,iostat=io_status)
      call assert( io_status==success, "initialize_input_parameters: "//pariz_namelist_file//" read." )
    end block

  end procedure

end submodule pariz_namelist_implementation
