module dBarrays
  use numeric_defs
  implicit none

  private
  protected :: DB, DDB, D3BX, D3BY, D3BZ, D4BX, D4BY, D4BZ

  real(dbl), dimension(3,3) :: DB
  real(dbl), dimension(3,3,3) :: DDB
  real(dbl), dimension(3,3,3) :: D3BX, D3BY, D3BZ
  real(dbl), dimension(3,3,3,3) :: D4BX, D4BY, D4BZ

  include "C.DDBXYZ.H"
  include "C.D3B_2.H"
  include "C.D4B.H"

end module dBarrays

