module taylor
  use numeric_defs
  implicit none

  private
  public :: ncdim, ntord
  public :: PascalMatrix, PascalEntry, PascalStart
  public :: orderStart, orderEnd, multinomCoeff, GIndexTimes
  public :: init_TSIndexing
  public :: GiorgilliIndex, GiorgilliExpon, stepGiorgilliExpon
  public :: evalMonoms, evalMonomsU


  integer, parameter :: ncdim = 3  ! num coordinate dimensions
  integer, parameter :: ntord = 5  ! max Taylor order

  ! PascalMatrix(0:ncdim, 0:ntord+1)
  !   = [ 0  1  1  1  1  1  1 ]
  !     [ 0  1  2  3  4  5  6 ]
  !     [ 0  1  3  6 10 15 21 ]
  !     [ 0  1  4 10 20 35 56 ]
  ! Matrix constructed from Pascal's triangle.
  ! All entries are binomial coefficients:
  !   PascalMatrix(j, k) = binom(j + k - 1, k - 1)
  ! or
  !   binom(n, k) = PascalMatrix(n - k, k + 1)
  integer, parameter :: PascalMatrix(0:ncdim, 0:ntord+1) = reshape( &
    & [ 0,  1,  1,  1,  1,  1,  1,   &
    &   0,  1,  2,  3,  4,  5,  6,   &
    &   0,  1,  3,  6, 10, 15, 21,   &
    &   0,  1,  4, 10, 20, 35, 56 ], &
    & [ ncdim + 1, ntord + 2 ], order = [2, 1])

  ! Pascal's triangle flattened into a linear array
  integer, allocatable :: PascalEntry(:)
  !
  ! index at which each row of Pascal's triangle begins
  ! ==> binom(n, k) = PascalEntry(PascalStart(n) + k)
  integer, parameter :: PascalStart(0:ntord+1) = PascalMatrix(ncdim-1,:)

  ! orderStart(0:ntord)
  !   = [ 0 1 4 10 20 35 ]
  ! Initial Giorgilli index for monomials of a given order.
  integer, parameter :: orderStart(0:ntord) = PascalMatrix(ncdim,0:ntord)

  ! orderEnd(0:ntord)
  !   = [ 0 3 9 19 34 55 ]
  ! Final Giorgilli index for monomials of a given order.
  integer, parameter :: orderEnd(0:ntord) = PascalMatrix(ncdim,1:ntord+1) - 1

  ! multinomCoeff(0:Smon(ntord,ncdim))
  !   = [ 1  1 1 1  1 2 2 1 2 1  1 3 3 3 6 3 1 3 3 1  ... ]
  ! Multinomial coefficients: These are a generalisation of
  ! the well-known binomial coefficients. For exponent vector
  !   [ j1 j2 j3 ], with j1 + j2 + j3 = n,
  ! we define the multinomial coefficient as
  !   n! / (j1! j2! j3!).
  !integer, allocatable :: multinomCoeff(:)
  integer, parameter :: multinomCoeff(0:orderEnd(ntord)) = &
    & [ 1, &
    &   1,  1,  1, &
    &   1,  2,  2,  1,  2,  1, &
    &   1,  3,  3,  3,  6,  3,  1,  3,  3,  1, &
    &   1,  4,  4,  6, 12,  6,  4, 12, 12,  4,  1,  4,  6,  4,  1, &
    &   1,  5,  5, 10, 20, 10, 10, 30, 30, 10,  5, 20, 30, 20,  5, &
    &           1,  5, 10, 10,  5,  1 ]

  ! look-up table for the Giorgilli index of a product of monomials
  type(iarray2) :: GIndexTimes


contains

  ! initialise the indexing arrays
  subroutine init_TSIndexing
    implicit none
    integer :: i, j, k, os
    integer :: gi, gi1, mc, nr
    integer :: fact(0:ntord)
    integer :: jv(1:ncdim)
    integer :: jv1(1:ncdim)

    ! initialise PascalEntry
    allocate(PascalEntry(0:((ntord+3) * (ntord+2) / 2 - 1)))
    ! allocate and populate the enrtries of Pascal's triangle
    do j = 0, ntord+1
      os = PascalStart(j)
      PascalEntry(os) = 1
      do k = 1, j-1
        PascalEntry(os + k) = PascalEntry(os - k)                       &
          &                 + PascalEntry(os - k - 1)
      end do
      PascalEntry(os + j) = 1
    end do

    !! initialise array of multinomial coefficients
    !allocate(multinomCoeff(0:orderEnd(ntord)))
    !fact(0) = 1
    !do i = 1, ntord
    !  fact(i) = i * fact(i-1)
    !end do
    !do i = 0, ntord
    !  do gi = orderStart(i), orderEnd(i)
    !    call GiorgilliExpon(i, gi, jv)
    !    mc = fact(i)
    !    do k = 1, ncdim
    !      mc = mc / fact(jv(k))
    !    end do
    !    multinomCoeff(gi) = mc
    !  end do
    !end do
    do i = 0, ntord
      write(6, '(55(i4))') (multinomCoeff(gi), gi=orderStart(i),orderEnd(i))
    end do

    ! allocate and populate the look-up table GIndexTimes
    allocate(GIndexTimes%c(0:orderEnd(ntord)))
    do i = 0, ntord
      nr = orderEnd(ntord - i)
      do gi = orderStart(i), orderEnd(i)
        allocate(GIndexTimes%c(gi)%r(0:nr))
        call GiorgilliExpon(i, gi, jv)
        jv1(:) = 0
        GIndexTimes%c(gi)%r(0) = gi
        do gi1 = 1, nr
          call stepGiorgilliExpon(jv1)
          GIndexTimes%c(gi)%r(gi1) = GiorgilliIndex(jv + jv1)
        end do
      end do
    end do
  end subroutine init_TSIndexing


  ! return Giorgilli index for a given exponent vector
  function GiorgilliIndex(jv) result(gi)
    implicit none
    integer, intent(in) :: jv(:)
    integer :: gi

    integer :: i, k, jsum, nv

    nv = size(jv)
    i = nv
    jsum = jv(i)
    gi = PascalMatrix(1, jsum)
    do k = 2, nv
      i = i - 1
      jsum = jsum + jv(i)
      gi = gi + PascalMatrix(k, jsum)
    end do
  end function GiorgilliIndex


  ! set exponent vector jv corresponding to Giorgilli index gi
  subroutine GiorgilliExpon(pord, gi, jv)
    implicit none
    integer, intent(in) :: pord   ! maximum order
    integer, intent(in) :: gi     ! Giorgilli index
    integer, intent(out) :: jv(:) ! resultant exponent vector (1:nv)

    integer :: gr, nv
    integer :: iord, iord1, iv, k

    nv = size(jv)
    gr = gi
    iord = pord
    iv = nv
    do while(PascalMatrix(iv, iord) .gt. gi)
      iord = iord - 1
    end do

    do k = 1, nv - 1
      iord1 = iord;
      gr = gr - PascalMatrix(iv, iord)
      iv = iv - 1
      do while(PascalMatrix(iv, iord) .gt. gr)
        iord = iord - 1
      end do
      jv(k) = iord1 - iord;
    end do
    jv(nv) = iord;
  end subroutine GiorgilliExpon


  ! advance exponent vector jv by one step in the Giorgilli sequence
  subroutine stepGiorgilliExpon(jv)
    implicit none
    integer, intent(inout) :: jv(:)
    integer :: nv, k, jn

    nv = size(jv)
    jn = jv(nv)
    jv(nv) = 0
    k = nv - 1
    do while(jv(k) .eq. 0 .and. k .gt. 0)
      k = k - 1
    end do
    if (k .ne. 0) jv(k) = jv(k) - 1
    jv(k + 1) = jn + 1
  end subroutine stepGiorgilliExpon


  ! Evaluate the sequence of monomials z^jv in Giorgilli order,
  ! returning mon <- z^0 .. z^pord  ([z] ~ nv, [mon] ~ nmon)
  ! with nmon = number of monomials of order pord in nv variables.
  subroutine evalMonoms(pord, z, mon)
    implicit none
    integer, intent(in) :: pord      ! maximum order
    real(dbl), intent(in) :: z(:)     ! linear monomials (1:nv)
    real(dbl), intent(out) :: mon(0:) ! resultant array of monomials (0:nmon)

    integer :: gi, nv
    integer :: k, m, ma, mai, maf

    nv = size(z)

    mon(0) = one;
    mon(1:nv) = z

    gi = nv + 1;
    do m = 2, pord
      maf = orderEnd(m - 1);
      do k = 1, nv
        mai = maf + 1 - PascalMatrix(nv - k, m)
        do ma = mai, maf
          mon(gi) = z(k) * mon(ma);
          gi = gi + 1
        end do
      end do
    end do
  end subroutine evalMonoms


  ! Evaluate the sequence of monomials in U, U', U'', etc., required
  ! to evaluate the derivatives B^(p) along a particle trajectory.
  ! These monomials include the multinomial coefficient as a factor,
  ! as well as additional numerical factors that depend on order p.
  ! This subroutine assumes it has been called previously for the
  ! orders 1 .. p-1. In addition, this subroutine assumes that the
  ! monomial array mon is either set to zero on first call, or never
  ! accessed beyond the entries relevant for computing B^(p).
  subroutine evalMonomsU(p, dU, mon)
    implicit none
    integer, intent(in) :: p          ! order of B-derivative
    real(dbl), intent(in) :: dU(:,0:)   ! dU(ivar, ideriv)
    real(dbl), intent(inout) :: mon(:) ! resultant array of U monomials (1:nmon)

    integer :: gi, k, ma, mai, maf

    select case (p)

    case (1) ! { U_i }
      mon(1:ncdim) = dU(:, 0)
      
    case (2) ! { U'_i, U_i U_j }
      ! populate the quadratic terms: U_i U_j
      gi = ncdim + 1; ! orderStart(2)
      maf = orderEnd(1);
      do k = 1, ncdim
        mai = maf + 1 - PascalMatrix(ncdim - k, 2)
        do ma = mai, maf
          mon(gi) = dU(k, 0) * mon(ma);
          gi = gi + 1
        end do
      end do
      mon(orderStart(2):orderEnd(2))                                    &
        & = mon(orderStart(2):orderEnd(2))                              &
        &   * multinomCoeff(orderStart(2):orderEnd(2))
      ! populate the linear terms: U'_i
      mon(1:ncdim) = dU(:, 1)
      
    case (3) ! { U''_i, 3 U'_i U_j, U_i U_j U_k }
      ! divide incoming quadratic terms by multinomial coefficients
      mon(orderStart(2):orderEnd(2))                                    &
        & = mon(orderStart(2):orderEnd(2))                              &
        &   / multinomCoeff(orderStart(2):orderEnd(2))
      ! populate the cubic terms: U_i U_j U_k
      gi = orderStart(3)
      maf = orderEnd(2);
      do k = 1, ncdim
        mai = maf + 1 - PascalMatrix(ncdim - k, 3)
        do ma = mai, maf
          mon(gi) = dU(k, 0) * mon(ma);
          gi = gi + 1
        end do
      end do
      mon(orderStart(3):orderEnd(3))                                    &
        & = mon(orderStart(3):orderEnd(3))                              &
        &   * multinomCoeff(orderStart(3):orderEnd(3))
      ! populate the quadratic terms: 3 U'_i U_j
      gi = orderStart(2)
      maf = orderEnd(1);
      do k = 1, ncdim
        mai = maf + 1 - PascalMatrix(ncdim - k, 3)
        do ma = mai, maf
          mon(gi) = dU(k, 0) * mon(ma);
          gi = gi + 1
        end do
      end do
      mon(5:6) = mon(5:6) + dU(1,1) * dU(2:3,0)
      mon(8) = mon(8) + dU(2,1) * dU(3,0)
      mon(orderStart(2):orderEnd(2))                                    &
        & = mon(orderStart(2):orderEnd(2)) * 3
      ! populate the linear terms: U''_i
      mon(1:ncdim) = dU(:, 2)
      
    case (4)
      ! { U'''_i, 4 U''_i U_j + 3 U'_i U'_j, 6 U'_i U_j U_k, U_i U_j U_k U_l }
      ! divide incoming cubic terms by multinomial coefficients
      mon(orderStart(3):orderEnd(3))                                    &
        & = mon(orderStart(3):orderEnd(3))                              &
        &   / multinomCoeff(orderStart(3):orderEnd(3))
      ! populate the quartic terms: U_i U_j U_k U_l
      gi = orderStart(4)
      maf = orderEnd(3);
      do k = 1, ncdim
        mai = maf + 1 - PascalMatrix(ncdim - k, 4)
        do ma = mai, maf
          mon(gi) = dU(k, 0) * mon(ma);
          gi = gi + 1
        end do
      end do
      mon(orderStart(4):orderEnd(4))                                    &
        & = mon(orderStart(4):orderEnd(4))                              &
        &   * multinomCoeff(orderStart(4):orderEnd(4))
      ! populate the cubic terms: 6 U'_i U_j U_k
      gi = orderStart(3)
      maf = orderEnd(2);
      do k = 1, ncdim
        mai = maf + 1 - PascalMatrix(ncdim - k, 3)
        do ma = mai, maf
          mon(gi) = dU(k, 0) * mon(ma);
          gi = gi + 1
        end do
      end do
      mon(11:12) = mon(11:12) + dU(2:3,0) * mon(4)
      mon(13:14) = mon(13:14) + dU(2:3,0) * mon(5)
      mon(14:15) = mon(14:15) + dU(2:3,0) * mon(6)
      mon(17:18) = mon(17:18) + dU(3,0) * mon(7:8)
      mon(orderStart(3):orderEnd(3))                                    &
        & = mon(orderStart(3):orderEnd(3)) * 2
      ! populate the quadratic terms: 4 U''_i U_j + 3 U'_i U'_j
      mon(4) = 4 *  dU(1,0) * dU(1,2)                      + 3 * dU(1,1) * dU(1,1)
      mon(5) = 4 * (dU(1,0) * dU(2,2) + dU(1,2) * dU(2,0)) + 6 * dU(1,1) * dU(2,1)
      mon(6) = 4 * (dU(1,0) * dU(3,2) + dU(1,2) * dU(3,0)) + 6 * dU(1,1) * dU(3,1)
      mon(7) = 4 *  dU(2,0) * dU(2,2)                      + 3 * dU(2,1) * dU(2,1)
      mon(8) = 4 * (dU(2,0) * dU(3,2) + dU(2,2) * dU(3,0)) + 6 * dU(2,1) * dU(3,1)
      mon(9) = 4 *  dU(3,0) * dU(3,2)                      + 3 * dU(3,1) * dU(3,1)
      ! populate the linear terms: U'''_i
      mon(1:ncdim) = dU(:, 3)
      
    case default
      ! emit fatal error
      write(*,*) ""
      write(*,*) "Fatal error!"
      write(*,*) "Subroutine evalMonomsU(p, dU, mon) in module taylor"
      write(*,'(a,i2,a)') "  called with argument p = ", p, "."
      write(*,*) "Exiting now!"
      error stop 1

    end select

  end subroutine evalMonomsU

end module taylor

