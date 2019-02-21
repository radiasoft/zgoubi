module particle
  use assertions_interface, only : assert, assertions
  use numeric_defs
  use taylor
  implicit none

  private
  public :: derivU                       !! Data
  public :: evalDU, DBarrays2dnB, init_particle_module !! Methods

  ! derivitives along the particle trajectory
  real(dbl) :: derivU(1:ncdim, 0:ntord)   ! normalized velocity u^(k)
  real(dbl) :: derivB(1:ncdim, 0:ntord-1) ! magnetic field      B^(k)

  ! spatial derivatives of the magnetic field B
  ! dnB_{x,y,z}, in Giorgilli order
  real(dbl), allocatable :: dnB(:, :)

  ! array of monomials involving U, U', etc.
  real(dbl), allocatable :: monomsU(:)

contains

  ! initialise memory for this module
  subroutine init_particle_module
    ! allocate dnB
    if(.not. allocated(dnB)) then
      allocate(dnB(ncdim, orderEnd(ntord-1)), source = zero_r)
    end if
    ! allocate monomsU
    if(.not. allocated(monomsU)) then
      allocate(monomsU(orderEnd(ntord-1)), source = zero_r)
    end if
  end subroutine init_particle_module


  ! compute d^(n+1) U / ds^(n+1)
  !   = sum_(k = 0 .. n) binom(n,k) * u^(n-k) x B^(k)
  subroutine evalDUn(np1, maxDB)
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
  end subroutine evalDUn


  ! compute d^(n+1) U / ds^(n+1) for n in 0:ntord
  ! NB: array dnB must be populated
  subroutine evalDU(U, B, maxDB)
    implicit none
    real(dbl), intent(in) :: U    ! normalized velocity U
    real(dbl), intent(in) :: B    ! normalized magnetic field B
    integer, intent(in) :: maxDB  ! maximuim order of B derivatives
                                  ! (e.g., uniform B => maxDB = 0)

    integer :: k

    ! Requires
    if (assertions) then
      call assert(lbound(derivB,2)==0 .and. ubound(derivB,ntord-1), "evalDU: received expected derivB shape")
    end if

    derivU = zero_r
    derivU(:,0) = U
    derivB(:,0) = B

    do k = 1, ntord-1
      call evalDUn(k,maxDB)
      call evalMonomsU(k, derivU, monomsU)
      derivB(:,k) = matmul( dnB(:, 1:orderEnd(k)), monomsU(1:orderEnd(k)) )
    end do
    call evalDUn(ntord, maxDB)

  end subroutine evalDU

  ! populate the one-dimensional field derivative array dnB
  ! from the various multi-dimensional DB arrays
  subroutine DBarrays2dnB(maxDB)
    use dBarrays
    implicit none
    integer, intent(in) :: maxDB  ! maximuim order of B derivatives
                                  ! (e.g., uniform B => maxDB = 0)

    integer :: i

    ! dB/dXi
    if(maxDB < 1) return
    dnB(:, 1:3) = DB(:,:)

    ! d^2 B / dXi.dXj
    if(maxDB < 2) return
    dnB(:, 4) = DDB(:, 1,1)
    dnB(:, 5) = DDB(:, 1,2)
    dnB(:, 6) = DDB(:, 1,3)
    dnB(:, 7) = DDB(:, 2,2)
    dnB(:, 8) = DDB(:, 2,3)
    dnB(:, 9) = DDB(:, 3,3)

    ! d^3 B / dXi.dXj.dXk
    if(maxDB < 3) return
    dnB(1, 10) = D3BX(1,1,1)
    dnB(2, 10) = D3BY(1,1,1)
    dnB(3, 10) = D3BZ(1,1,1)
    dnB(1, 11) = D3BX(1,1,2)
    dnB(2, 11) = D3BY(1,1,2)
    dnB(3, 11) = D3BZ(1,1,2)
    dnB(1, 12) = D3BX(1,1,3)
    dnB(2, 12) = D3BY(1,1,3)
    dnB(3, 12) = D3BZ(1,1,3)
    dnB(1, 13) = D3BX(1,2,2)
    dnB(2, 13) = D3BY(1,2,2)
    dnB(3, 13) = D3BZ(1,2,2)
    dnB(1, 14) = D3BX(1,2,3)
    dnB(2, 14) = D3BY(1,2,3)
    dnB(3, 14) = D3BZ(1,2,3)
    dnB(1, 15) = D3BX(1,3,3)
    dnB(2, 15) = D3BY(1,3,3)
    dnB(3, 15) = D3BZ(1,3,3)
    dnB(1, 16) = D3BX(2,2,2)
    dnB(2, 16) = D3BY(2,2,2)
    dnB(3, 16) = D3BZ(2,2,2)
    dnB(1, 17) = D3BX(2,2,3)
    dnB(2, 17) = D3BY(2,2,3)
    dnB(3, 17) = D3BZ(2,2,3)
    dnB(1, 18) = D3BX(2,3,3)
    dnB(2, 18) = D3BY(2,3,3)
    dnB(3, 18) = D3BZ(2,3,3)
    dnB(1, 19) = D3BX(3,3,3)
    dnB(2, 19) = D3BY(3,3,3)
    dnB(3, 19) = D3BZ(3,3,3)

    ! d^4 B / dXi.dXj.dXk.dXl
    if(maxDB < 4) return
    dnB(1, 20) = D4BX(1,1,1,1)
    dnB(2, 20) = D4BY(1,1,1,1)
    dnB(3, 20) = D4BZ(1,1,1,1)
    dnB(1, 21) = D4BX(1,1,1,2)
    dnB(2, 21) = D4BY(1,1,1,2)
    dnB(3, 21) = D4BZ(1,1,1,2)
    dnB(1, 22) = D4BX(1,1,1,3)
    dnB(2, 22) = D4BY(1,1,1,3)
    dnB(3, 22) = D4BZ(1,1,1,3)
    dnB(1, 23) = D4BX(1,1,2,2)
    dnB(2, 23) = D4BY(1,1,2,2)
    dnB(3, 23) = D4BZ(1,1,2,2)
    dnB(1, 24) = D4BX(1,1,2,3)
    dnB(2, 24) = D4BY(1,1,2,3)
    dnB(3, 24) = D4BZ(1,1,2,3)
    dnB(1, 25) = D4BX(1,1,3,3)
    dnB(2, 25) = D4BY(1,1,3,3)
    dnB(3, 25) = D4BZ(1,1,3,3)
    dnB(1, 26) = D4BX(1,2,2,2)
    dnB(2, 26) = D4BY(1,2,2,2)
    dnB(3, 26) = D4BZ(1,2,2,2)
    dnB(1, 27) = D4BX(1,2,2,3)
    dnB(2, 27) = D4BY(1,2,2,3)
    dnB(3, 27) = D4BZ(1,2,2,3)
    dnB(1, 28) = D4BX(1,2,3,3)
    dnB(2, 28) = D4BY(1,2,3,3)
    dnB(3, 28) = D4BZ(1,2,3,3)
    dnB(1, 29) = D4BX(1,3,3,3)
    dnB(2, 29) = D4BY(1,3,3,3)
    dnB(3, 29) = D4BZ(1,3,3,3)
    dnB(1, 30) = D4BX(2,2,2,2)
    dnB(2, 30) = D4BY(2,2,2,2)
    dnB(3, 30) = D4BZ(2,2,2,2)
    dnB(1, 31) = D4BX(2,2,2,3)
    dnB(2, 31) = D4BY(2,2,2,3)
    dnB(3, 31) = D4BZ(2,2,2,3)
    dnB(1, 32) = D4BX(2,2,3,3)
    dnB(2, 32) = D4BY(2,2,3,3)
    dnB(3, 32) = D4BZ(2,2,3,3)
    dnB(1, 33) = D4BX(2,3,3,3)
    dnB(2, 33) = D4BY(2,3,3,3)
    dnB(3, 33) = D4BZ(2,3,3,3)
    dnB(1, 34) = D4BX(3,3,3,3)
    dnB(2, 34) = D4BY(3,3,3,3)
    dnB(3, 34) = D4BZ(3,3,3,3)
  end subroutine DBarrays2dnB


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

