submodule(particle) particle_impl
  implicit none

contains

  module procedure evalDU
    use taylor, only : PascalStart, PascalEntry, ncdim
    use numeric_defs, only : dbl, zero_r
    integer :: k, km, n, psn
    real(dbl) :: uxb(1:ncdim)

    n = np1 - 1
    psn = PascalStart(n)
    km = min(n, maxDB)
    derivU(:,np1) = zero_r
    do k = 0, km
      uxb(:) = derivU(:,n-k) .cross. derivB(:,k)
      derivU(:,np1) = derivU(:,np1) + PascalEntry(psn + k) * uxb(:)
    end do
  end procedure evalDU


  module procedure derivU_column
    use taylor, only : PascalStart, PascalEntry, ncdim
    use numeric_defs, only : dbl, zero_r
    integer :: k

    associate( n => np1 - 1 )
    associate( psn => PascalStart(n) )
    associate( km => min(n, maxDB) )
      dUn = sum([(PascalEntry(psn + k) * (derivU(:,n-k) .cross. derivB(:,k)), k=0,km)])
    end associate; end associate; end associate

  end procedure derivU_column


  module procedure cross
    axb(1) = a(2) * b(3) - a(3) * b(2)
    axb(2) = a(3) * b(1) - a(1) * b(3)
    axb(3) = a(1) * b(2) - a(2) * b(1)
  end procedure cross

end submodule particle_impl
