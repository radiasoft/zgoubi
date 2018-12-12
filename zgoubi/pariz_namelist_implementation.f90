submodule(pariz_namelist_interface) pariz_namelist_implementation
  use assertions_interface, only : assert
  implicit none

contains

  module procedure initialize_input_parameters

    namelist/pariz/ ID, IZ, MMAP, MXX, MXY

    if ( all( [ID,IZ,MMAP,MXX,MXY]/=unread) ) return  !! Return if the namelist has already been read

    block 
      integer, parameter :: success=0
      integer file_unit, io_status
      logical pariz_nml_open

      inquire(file=pariz_namelist_file, number=file_unit, opened=pariz_nml_open)
      call assert( .not. pariz_nml_open, "read_from_file: file not already open" )
        !! If the above assertion is removed, then add code the following form:
        !! if(parize_nml_open) then; rewind(...); else; open(...)
    
      open(newunit=file_unit,file=pariz_namelist_file,iostat=io_status,status='old')
      call assert( io_status==success, "input_initialization_parameters: "//pariz_namelist_file//" opened." )
  
      read(file_unit,nml=pariz,iostat=io_status)
      call assert( io_status==success, "input_initialization_parameters: "//pariz_namelist_file//" read." )
    end block

  end procedure

  module procedure all_namelist_values_read
    call assert( all( [ID,IZ,MMAP,MXX,MXY]/=unread ), "input_initialization_parameters: all namelist values read" )
  end procedure 

end submodule pariz_namelist_implementation
