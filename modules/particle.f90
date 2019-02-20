module particle
  use numeric_defs
  use taylor
  implicit none

  private
  public :: derivU, derivB, evalDU

  ! derivitives along the particle trajectory
  real(dbl) :: derivU(1:ncdim, 0:ntord)   ! normalized velocity u^(k)
  real(dbl) :: derivB(1:ncdim, 0:ntord-1) ! magnetic field      B^(k)

contains

  ! compute d^(n+1) U / ds^(n+1)
  !   = sum_(k = 0 .. n) binom(n,k) * u^(n-k) x B^(k)
  subroutine evalDU(np1, maxDB)
    implicit none
    integer, intent(in) :: np1    ! derivative order to evaluate (n + 1)
    integer, intent(in) :: maxDB  ! maximuim order of B derivatives
                                  ! (e.g., uniform B => maxDB = 0)

    integer :: k, km, n, psn
    real(dbl) :: uxb(1:ncdim)

    n = np1 - 1
    psn = PascalStart(n)
    km = min(n, maxDB)
    derivU(:,np1) = zero_r
    do k = 0, km
      uxb(:) = cross(derivU(:,n-k), derivB(:,k))
      derivU(:,np1) = derivU(:,np1) + PascalEntry(psn + k) * uxb(:)
    end do
  end subroutine evalDU


  ! return the standard 3-vector cross-product a x b
  function cross(a, b)  result(axb)
    implicit none
    real(dbl), intent(in) :: a(3)  ! left-hand argument
    real(dbl), intent(in) :: b(3)  ! right-hand argument
    real(dbl) :: axb(3)
    axb(1) = a(2) * b(3) - a(3) * b(2)
    axb(2) = a(3) * b(1) - a(1) * b(3)
    axb(3) = a(1) * b(2) - a(2) * b(1)
  end function cross


end module particle

