!> \file constants.F90
!! \BRIEF
!> Module with contants subroutine - computes carbonate system constants
!! from S,T,P
module mocsy_physical_constants

  use mocsy_singledouble, only: rx, r8, wp

  implicit none; private

  ! physical constants 
  real(r8), parameter, public :: ideal_gas_constant_jkmol = 8.314472_r8 ! [J/(mol*K)]
  real(r8), parameter, public :: R_jkmol_scaled_by_10 = 83.14472_r8
  real(r8), parameter, public :: ideal_gas_constant_codata = 82.05736_r8 ! [cm^3*atm/(K*mol)]
  real(r8), parameter, public :: co2_partial_molar_volume = 32.3_r8  ! [cm3/mol]
  real(r8), parameter, public :: zero_c_in_kelvin = 273.15_r8
  real(r8), parameter, public :: bar_to_atm = 1.01325_r8
  real(r8), parameter, public :: fluoride_atomic_mass = 18.9984_r8! [F-]
  real(r8), parameter, public :: sulfate_atomic_mass = 96.062_r8! [SO4-]
  real(r8), parameter, public :: boron_atomic_mass = 10.811_r8 ! [B]
  real(r8), parameter, public :: knudsen_chlorinity_constant = 1.80655_r8 ! based on Knudsen's formula
real(r8), parameter, public :: pi = 3.141592654_r8 ! we should make this better


end module mocsy_physical_constants
