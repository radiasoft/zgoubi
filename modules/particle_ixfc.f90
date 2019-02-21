module particle

  use numeric_defs, only : dbl
  use taylor, only : ncdim, ntord
  implicit none

  private
  public :: derivU, derivB, evalDU
  public :: operator(.cross.)

  real(dbl) :: derivU(1:ncdim, 0:ntord)   !! normalized velocity u^(k)
   !! normalized velocity derivatives, u^(k), along the particle trajectory
  real(dbl) :: derivB(1:ncdim, 0:ntord-1) 
   !! magnetic field derivatives, B^(k), along the particle trajectory

  interface operator(.cross.)
    !! cross-product operator
    module procedure cross
  end interface

  interface 

    module subroutine evalDU(np1, maxDB)
      !! compute d^(n+1) U / ds^(n+1)
      !!   = sum_(k = 0 .. n) binom(n,k) * u^(n-k) x B^(k)
      !! (e.g., uniform B => maxDB = 0)
      implicit none
      integer, intent(in) :: np1    !! derivative order to evaluate (n + 1)
      integer, intent(in) :: maxDB  !! maximuim order of B derivatives
    end subroutine evalDU

    pure module function derivU_column(np1, maxDB) result(dUn)
      !! compute d^(n+1) U / ds^(n+1)
      !!   = sum_(k = 0 .. n) binom(n,k) * u^(n-k) x B^(k)
      !! (e.g., uniform B => maxDB = 0)
      implicit none
      integer, intent(in) :: np1    !! derivative order to evaluate (n + 1)
      integer, intent(in) :: maxDB  !! maximuim order of B derivatives
      real(dbl) :: dUn(1:ncdim)   !! normalized velocity u^(k)
    end function derivU_column

    pure module function cross(a, b)  result(axb)
      !! return the standard 3-vector cross-product a x b
      implicit none
      real(dbl), intent(in) :: a(3)  !! left-hand argument
      real(dbl), intent(in) :: b(3)  !! right-hand argument
      real(dbl) :: axb(3)            !! cross-product
    end function cross

  end interface

end module particle

