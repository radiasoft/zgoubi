module mathphys_consts
  use numeric_defs, only : dbl, half, one, two, four
  implicit none

  private
   !! Nothing in this module is used yet.
   !! Let's explicitly add names to a public list as necessary.

  !
  ! mathematical constants
  !
  !real(dbl), parameter :: pi = 3.141592653589793238462643383279502e0_dbl
  real(dbl), parameter :: pi = four * atan(one)
  real(dbl), parameter :: twopi = two * pi
  real(dbl), parameter :: fourpi = four * pi
  real(dbl), parameter :: halfpi = half * pi
  real(dbl), parameter :: invtwopi = one / twopi
  real(dbl), parameter :: invfourpi = one / fourpi
  real(dbl), parameter :: euler_e = exp(one)
  real(dbl), parameter :: rad_to_deg = 180.0_dbl / pi
  real(dbl), parameter :: deg_to_rad = pi / 180.0_dbl

  !
  ! physical constants: 2014 CODATA recommended values
  ! http://physics.nist.gov/cuu/index.html
  ! The uncertainties (given in parentheses) indicate
  ! the expected error in the last digits.
  !
  ! universal
  !
  ! speed of light in vacuum / m.s^-1  (exact)
  real(dbl), parameter :: c_light = 2.99792458e8_dbl
  ! magnetic constant (= permeability of free space) / H.m^-1  (exact)
  real(dbl), parameter :: mu_0 = 12.5663706143591729539e-7_dbl
  ! {magnetic constant} / (4.pi) / H.m^-1  (exact)  [H.m^-1 == N.A^-2]
  real(dbl), parameter :: mu_0_ovr_4pi = 1.e-7_dbl
  ! electric constant (= permittivity of free space) / F.m^-1  (exact)
  real(dbl), parameter :: epsilon_0 = 8.85418781762038985054e-12_dbl
  ! 4.pi.{electric constant} / F.m^-1  (exact)
  real(dbl), parameter :: four_pi_epsilon_0 = 1.11265005605361843217e-10_dbl
  ! Planck constant / J.s  (81)
  real(dbl), parameter :: planck_h = 6.626070040e-34_dbl
  ! {Planck constant} / (2.pi) / J.s  (13)
  real(dbl), parameter :: planck_hbar = 1.054571800e-34_dbl
  ! Planck constant / eV.s  (25)
  real(dbl), parameter :: planck_h_evs = 4.135667662e-15_dbl
  ! {Planck constant} / (2.pi) /  eV.s  (40)
  real(dbl), parameter :: planck_hbar_evs = 6.582119514e-16_dbl
  ! {Planck constant}.c / (2.pi) /  eV.m  (12)
  real(dbl), parameter :: hbar_clight_evm = 197.3269788e-9_dbl
  !
  ! electromagnetic
  !
  ! elementary charge / C  (98)
  real(dbl), parameter :: elem_charge = 1.6021766208e-19_dbl
  ! Bohr magneton / J.T^-1  (57)
  real(dbl), parameter :: bohr_magneton = 927.4009994e-26_dbl
  ! Bohr magneton / eV.T^-1  (26)
  real(dbl), parameter :: bohr_magneton_ev = 5.7883818012e-5_dbl
  ! nuclear magneton / J.T^-1  (31)
  real(dbl), parameter :: nuclear_magneton = 5.050783699e-27_dbl
  ! nuclear magneton / eV.T^-1  (15)
  real(dbl), parameter :: nuclear_magneton_ev = 3.1524512550e-8_dbl
  !
  ! atomic and nuclear
  !
  ! electron mass / kg  (11)
  real(dbl), parameter :: electron_mass = 9.10938356e-31_dbl
  ! electron rest energy / eV  (31)
  real(dbl), parameter :: electron_restE = 0.5109989461e6_dbl
  ! muon mass / kg  (48)
  real(dbl), parameter :: muon_mass = 1.883531594e-28_dbl
  ! muon rest energy / eV  (24)
  real(dbl), parameter :: muon_restE = 105.6583745e6_dbl
  ! proton mass / kg  (21)
  real(dbl), parameter :: proton_mass = 1.672621898e-27_dbl
  ! proton rest energy / eV  (58)
  real(dbl), parameter :: proton_restE = 938.2720813e6_dbl
  ! deuteron mass / kg  (41)
  real(dbl), parameter :: deuteron_mass = 3.343583719e-27_dbl
  ! deuteron rest energy / eV  (12)
  real(dbl), parameter :: deuteron_restE = 1875.612928e6_dbl
  ! electron gyromagnetic anomaly  (26)
  real(dbl), parameter :: electron_G = 1.15965218091e-3_dbl
  ! muon gyromagnetic anomaly  (63)
  real(dbl), parameter :: muon_G = 1.16592089e-3_dbl
  ! proton gyromagnetic anomaly  (85)
  real(dbl), parameter :: proton_G = 1.7928473510_dbl
  ! deuteron gyromagnetic anomaly  (21)
  real(dbl), parameter :: deuteron_G = -0.142987272_dbl
  ! fine-structure constant  (17)
  real(dbl), parameter :: fine_structure_a = 7.2973525664e-3_dbl
  ! fine-structure constant  (31)
  real(dbl), parameter :: inv_fine_structure_a = 137.035999139_dbl
  ! Bohr radius / m  (12)
  real(dbl), parameter :: bohr_radius = 0.52917721067e-10_dbl
  ! classical electron radius / m  (19)
  real(dbl), parameter :: electron_radius = 2.8179403227e-15_dbl
  !
  ! physico-chemical
  !
  ! Boltzmann constant / J.K^-1  (79)
  real(dbl), parameter :: boltzmann_k = 1.38064852e-23_dbl
  ! Boltzmann constant / eV.K^-1  (50)
  real(dbl), parameter :: boltzmann_k_ev = 8.6173303e-5_dbl

end module mathphys_consts

