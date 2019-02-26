submodule(particle) particle_impl
  use assertions_interface, only : assert, assertions
  use numeric_defs, only : zero_r
  use taylor, only : orderEnd
  implicit none

  real(dbl) :: derivB(1:ncdim, 0:ntord-1)
   !! magnetic field derivatives, B^(k), along the particle trajectory

  real(dbl), allocatable :: dnB(:, :)
    !! spatial derivatives of the magnetic field B
    !! dnB_{x,y,z}, in Giorgilli order

  real(dbl), allocatable :: monomsU(:)
    !! array of monomials involving U, U', etc.

contains

  module procedure init_particle_module
    !! initialise memory for this module
    !! allocate dnB
    if(.not. allocated(dnB)) allocate(dnB(ncdim, orderEnd(ntord-1)), source = zero_r)
    !! allocate monomsU
    if(.not. allocated(monomsU)) allocate(monomsU(orderEnd(ntord-1)), source = zero_r)
  end procedure init_particle_module

  module procedure evalDUn
    use taylor, only : PascalStart, PascalEntry, ncdim
    use numeric_defs, only : dbl, zero_r
    integer :: k, km, n, psn
    real(dbl) :: uxb(1:ncdim)

    n = np1 - 1
    psn = PascalStart(n)
    km = min(n, maxDB)

    if (assertions) then
      call assert(all([lbound(derivU,2) <= np1 , np1  <= ubound(derivU,2)]), "evalDUn: derivU(:,np1) in bounds" )
      call assert(all([lbound(derivU,2) <= n   , n    <= ubound(derivU,2)]), "evalDUn: derivU(:,n) in bounds"   )
      call assert(all([lbound(derivU,2) <= n-km, n-km <= ubound(derivU,2)]), "evalDUn: derivU(:,n-km) in bounds")
      call assert(all([lbound(derivB,2) <= 0   , 0    <= ubound(derivB,2)]), "evalDUn: derivB(:,0) in bounds"   )
      call assert(all([lbound(derivB,2) <= km  , km   <= ubound(derivB,2)]), "evalDUn: derivB(:,km) in bounds"  )
      call assert(all([lbound(PascalEntry,1) <= psn,    psn    <= ubound(PascalEntry,1)]), "evalDUn: PascalEntry(psn + k) in bounds")
      call assert(all([lbound(PascalEntry,1) <= psn+km, psn+km <= ubound(PascalEntry,1)]), "evalDUn: PascalEntry(psn + k) in bounds")
      call assert(all([size(derivU,1),size(derivB,1)] == size(uxb)), "evalDUn: conformable cross-product operands & result")
    end if

    derivU(:,np1) = zero_r

    do k = 0, km
      uxb(:) = derivU(:,n-k) .cross. derivB(:,k)
      derivU(:,np1) = derivU(:,np1) + PascalEntry(psn + k) * uxb(:)
    end do
  end procedure evalDUn


  module procedure evalDU
    use numeric_defs, only : zero_r
    use taylor, only : orderEnd, evalMonomsU
    integer :: k

    ! Requires
    if (assertions) then
      call assert(lbound(derivB,2)==0 .and. ubound(derivB,2)==ntord-1, "evalDU: received expected derivB shape")
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
  end procedure evalDU


  module procedure derivU_column
    use taylor, only : PascalStart, PascalEntry
    integer :: k

    associate( n => np1 - 1 )
    associate( psn => PascalStart(n) )
    associate( km => min(n, maxDB) )
      dUn = sum([(PascalEntry(psn + k) * (derivU(:,n-k) .cross. derivB(:,k)), k=0,km)])
    end associate; end associate; end associate
  end procedure derivU_column


  module procedure DBarrays2dnB
    use dBarrays, only : DB, DDB, D3BX, D3BY, D3BZ, D4BX, D4BY, D4BZ

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
  end procedure DBarrays2dnB


  module procedure cross
    axb(1) = a(2) * b(3) - a(3) * b(2)
    axb(2) = a(3) * b(1) - a(1) * b(3)
    axb(3) = a(1) * b(2) - a(2) * b(1)
  end procedure cross

end submodule particle_impl
