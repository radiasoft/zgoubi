module dBarrays
  use numeric_defs
  implicit none

  private
  public :: DB, DDB, D3BX, D3BY, D3BZ, D4BX, D4BY, D4BZ

  real(dbl) :: DB
  real(dbl) :: DDB
  real(dbl) :: D3BX, D3BY, D3BZ
  real(dbl) :: D4BX, D4BY, D4BZ

  include "C.DDBXYZ.H"
  include "C.D3B_2.H"
  include "C.D4B.H"

end module dBarrays

