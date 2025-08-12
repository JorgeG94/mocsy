!> \file constants.F90
!! \BRIEF
!> Module with contants subroutine - computes carbonate system constants
!! from S,T,P
module mocsy_constants

   use mocsy_singledouble, only : rx, r8, wp
   use mocsy_p80, only : p80
   use mocsy_sw_temp, only : sw_temp, sw_temp_DNAD
   use mocsy_sw_ptmp, only : sw_ptmp, sw_ptmp_DNAD
   use Dual_Num_Auto_Diff

   implicit none ; private

   public constants

   real(r8), parameter, public :: ideal_gas_constant_jkmol = 8.314472_r8 ! [J/(mol*K)]
   real(r8), parameter, public :: R_jkmol_scaled_by_10 = 83.14472_r8
   real(r8), parameter, public :: ideal_gas_constant_codata = 82.05736_r8 ! [cm^3*atm/(K*mol)]
   real(r8), parameter, public :: co2_partial_molar_volume = 32.3_r8  ! [cm3/mol]
   real(r8), parameter, public :: zero_c_in_kelvin = 273.15_r8
   real(r8), parameter, public :: bar_to_atm =  1.01325_r8
   real(r8), parameter, public :: fluoride_atomic_mass = 18.9984_r8! [F-]
   real(r8), parameter, public :: sulfate_atomic_mass = 96.062_r8! [SO4-]
   real(r8), parameter, public :: boron_atomic_mass = 10.811_r8 ! [B]
   real(r8), parameter, public :: knudsen_chlorinity_constant = 1.80655_r8 ! based on Knudsen's formula

   ! CONSTANTS
   ! =========
   ! Constants in formulation for Pressure effect on K's (Millero, 95)
   ! with corrected coefficients for Kb, Kw, Ksi, etc.

   ! index: 1) K1 , 2) K2, 3) Kb, 4) Kw, 5) Ks, 6) Kf, 7) Kspc, 8) Kspa,
   !            9) K1P, 10) K2P, 11) K3P, 12) Ksi


   real(r8), parameter :: a0(12) = [-25.5_r8, -15.82_r8, -29.48_r8, -20.02_r8, &
      -18.03_r8,  -9.78_r8, -48.76_r8, -45.96_r8, &
      -14.51_r8, -23.12_r8, -26.57_r8, -29.48_r8]

   real(r8), parameter :: a1(12) = [0.1271_r8, -0.0219_r8, 0.1622_r8, 0.1119_r8, &
      0.0466_r8, -0.0090_r8, 0.5304_r8, 0.5304_r8, &
      0.1211_r8, 0.1758_r8, 0.2020_r8, 0.1622_r8]

   real(r8), parameter :: a2(12) = [0.0_r8, 0.0_r8, -2.608e-3_r8, -1.409e-3_r8, &
      0.316e-3_r8, -0.942e-3_r8, 0.0_r8, 0.0_r8, &
      -0.321e-3_r8, -2.647e-3_r8, -3.042e-3_r8, -2.6080e-3_r8]

   real(r8), parameter :: b0(12) = [-3.08e-3_r8, 1.13e-3_r8, -2.84e-3_r8, -5.13e-3_r8, &
      -4.53e-3_r8, -3.91e-3_r8, -11.76e-3_r8, -11.76e-3_r8, &
      -2.67e-3_r8, -5.15e-3_r8, -4.08e-3_r8, -2.84e-3_r8]

   real(r8), parameter :: b1(12) = [0.0877e-3_r8, -0.1475e-3_r8, 0.0_r8, 0.0794e-3_r8, &
      0.09e-3_r8, 0.054e-3_r8, 0.3692e-3_r8, 0.3692e-3_r8, &
      0.0427e-3_r8, 0.09e-3_r8, 0.0714e-3_r8, 0.0_r8]

   real(r8), parameter :: b2(12) = 0.0_r8  ! All elements set to 0.0



contains
!> Compute thermodynamic constants
!! FROM temperature, salinity, and pressure (1D arrays)
   subroutine constants(k0, k1, k2, Kb, Kw, Ks, Kf, Kspc, Kspa,  &
      K1p, K2p, K3p, Ksi,                      &
      St, Ft, Bt,                              &
      temp, sal, Patm,                         &
      depth, lat, number_of_records,                           &
      optT, optP, optB, optK1K2, optKf, optGAS)

      !   Purpose:
      !     Compute thermodynamic constants
      !     FROM: temperature, salinity, and pressure (1D arrays)

      !     INPUT variables:
      !     ================
      !     Patm    = atmospheric pressure [atm]
      !     depth   = depth [m]     (with optP='m', i.e., for a z-coordinate model vertical grid is depth, not pressure)
      !             = pressure [db] (with optP='db')
      !     lat     = latitude [degrees] (needed to convert depth to pressure, i.e., when optP='m')
      !             = dummy array (unused when optP='db')
      !     temp    = potential temperature [degrees C] (with optT='Tpot', i.e., models carry tempot, not temp)
      !             = in situ   temperature [degrees C] (with optT='Tinsitu', e.g., for data)
      !     sal     = salinity in [psu]
      !     ---------
      !     optT: choose in situ vs. potential temperature as input
      !     ---------
      !     NOTE: Carbonate chem calculations require IN-SITU temperature (not potential Temperature)
      !       -> 'Tpot' means input is pot. Temperature (in situ Temp "tempis" is computed)
      !       -> 'Tinsitu' means input is already in-situ Temperature, not pot. Temp ("tempis" not computed)
      !     ---------
      !     optP: choose depth (m) vs pressure (db) as input
      !     ---------
      !       -> 'm'  means "depth" input is in "m" (thus in situ Pressure "p" [db] is computed)
      !       -> 'db' means "depth" input is already in situ pressure [db], not m (p = depth)
      !     ---------
      !     optB:
      !     ---------
      !       -> 'u74' means use classic formulation of Uppström (1974) for total Boron
      !       -> 'l10' means use newer   formulation of Lee et al. (2010) for total Boron
      !     ---------
      !     optK1K2:
      !     ---------
      !       -> 'l'   means use Lueker et al. (2000) formulations for K1 & K2 (recommended by Dickson et al. 2007)
      !                **** BUT this should only be used when 2 < T < 35 and 19 < S < 43
      !       -> 'm10' means use Millero (2010) formulation for K1 & K2 (see Dickson et al., 2007)
      !                **** Valid for 0 < T < 50°C and 1 < S < 50 psu
      !       -> 'w14' means use Waters (2014) formulation for K1 & K2 (see Dickson et al., 2007)
      !                **** Valid for 0 < T < 50°C and 1 < S < 50 psu
      !     -----------
      !     optKf:
      !     ----------
      !       -> 'pf' means use Perez & Fraga (1987) formulation for Kf (recommended by Dickson et al., 2007)
      !               **** BUT Valid for  9 < T < 33°C and 10 < S < 40.
      !       -> 'dg' means use Dickson & Riley (1979) formulation for Kf (recommended by Dickson & Goyet, 1994)
      !     -----------
      !     optGAS: choose in situ vs. potential fCO2 and pCO2
      !     ---------
      !       PRESSURE corrections for K0 and the fugacity coefficient (Cf)
      !       -> 'Pzero'   = 'zero order' fCO2 and pCO2 (typical approach, which is flawed)
      !                      considers in situ T & only atm pressure (hydrostatic=0)
      !       -> 'Ppot'    = 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
      !                      considers potential T & only atm pressure (hydrostatic press = 0)
      !       -> 'Pinsitu' = 'in situ' fCO2 and pCO2 (accounts for huge effects of pressure)
      !                      considers in situ T & total pressure (atm + hydrostatic)
      !     ---------

      !     OUTPUT variables:
      !     =================
      !     K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi
      !     St, Ft, Bt


#if USE_PRECISION == 2
#   define SGLE(x)    (x)
#else
#   define SGLE(x)    REAL(x)
#endif
! Input variables
      !>     number of records
!f2py intent(hide) number_of_records
      integer, intent(in) :: number_of_records
      !> in <b>situ temperature</b> (when optT='Tinsitu', typical data)
      !! OR <b>potential temperature</b> (when optT='Tpot', typical models) [degree C]
      real(kind=rx), intent(in),    dimension(number_of_records) :: temp
      !> depth in <b>meters</b> (when optP='m') or <b>decibars</b> (when optP='db')
      real(kind=rx), intent(in),    dimension(number_of_records) :: depth
      !> latitude <b>[degrees north]</b>
      real(kind=rx), intent(in),    dimension(number_of_records) :: lat
      !> salinity <b>[psu]</b>
      real(kind=rx), intent(in), dimension(number_of_records) :: sal

      !> atmospheric pressure <b>[atm]</b>
      real(kind=rx), intent(in), dimension(number_of_records) :: Patm

      !> for temp input, choose \b 'Tinsitu' for in situ Temp or
      !! \b 'Tpot' for potential temperature (in situ Temp is computed, needed for models)
      character(7), intent(in) :: optT
      !> for depth input, choose \b "db" for decibars (in situ pressure) or \b "m" for meters (pressure is computed, needed for models)
      character(2), intent(in) :: optP
      !> for total boron, choose either \b 'u74' (Uppstrom, 1974) or \b 'l10' (Lee et al., 2010).
      !! The 'l10' formulation is based on 139 measurements (instead of 20),
      !! uses a more accurate method, and
      !! generally increases total boron in seawater by 4%
!f2py character*3 optional, intent(in) :: optB='l10'
      character(3), optional, intent(in) :: optB
      !> for Kf, choose either \b 'pf' (Perez & Fraga, 1987) or \b 'dg' (Dickson & Riley, 1979)
!f2py character*2 optional, intent(in) :: optKf='pf'
      character(2), optional, intent(in) :: optKf
      !> for K1,K2 choose either \b 'l' (Lueker et al., 2000) or \b 'm10' (Millero, 2010) or \b 'w14' (Waters et al., 2014)
!f2py character*3 optional, intent(in) :: optK1K2='l'
      character(3), optional, intent(in) :: optK1K2
      !> for K0,fugacity coefficient choose either \b 'Ppot' (no pressure correction) or \b 'Pinsitu' (with pressure correction)
      !! 'Ppot'    - for 'potential' fCO2 and pCO2 (water parcel brought adiabatically to the surface)
      !! 'Pinsitu' - for 'in situ' values of fCO2 and pCO2, accounting for pressure on K0 and Cf
      !! with 'Pinsitu' the fCO2 and pCO2 will be many times higher in the deep ocean
!f2py character*7 optional, intent(in) :: optGAS='Pinsitu'
      character(7), optional, intent(in) :: optGAS


! Ouput variables
      !> solubility of CO2 in seawater (Weiss, 1974), also known as K0
      real(kind=r8), intent(out), dimension(number_of_records) :: k0
      !> K1 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
      real(kind=r8), intent(out), dimension(number_of_records) :: k1
      !> K2 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
      real(kind=r8), intent(out), dimension(number_of_records) :: k2
      !> equilibrium constant for dissociation of boric acid
      real(kind=r8), intent(out), dimension(number_of_records) :: Kb
      !> equilibrium constant for the dissociation of water (Millero, 1995)
      real(kind=r8), intent(out), dimension(number_of_records) :: Kw
      !> equilibrium constant for the dissociation of bisulfate (Dickson, 1990)
      real(kind=r8), intent(out), dimension(number_of_records) :: Ks
      !> equilibrium constant for the dissociation of hydrogen fluoride
      !! either from Dickson and Riley (1979) or from Perez and Fraga (1987), depending on optKf
      real(kind=r8), intent(out), dimension(number_of_records) :: Kf
      !> solubility product for calcite (Mucci, 1983)
      real(kind=r8), intent(out), dimension(number_of_records) :: Kspc
      !> solubility product for aragonite (Mucci, 1983)
      real(kind=r8), intent(out), dimension(number_of_records) :: Kspa
      !> 1st dissociation constant for phosphoric acid (Millero, 1995)
      real(kind=r8), intent(out), dimension(number_of_records) :: K1p
      !> 2nd dissociation constant for phosphoric acid (Millero, 1995)
      real(kind=r8), intent(out), dimension(number_of_records) :: K2p
      !> 3rd dissociation constant for phosphoric acid (Millero, 1995)
      real(kind=r8), intent(out), dimension(number_of_records) :: K3p
      !> equilibrium constant for the dissociation of silicic acid (Millero, 1995)
      real(kind=r8), intent(out), dimension(number_of_records) :: Ksi
      !> total sulfate (Morris & Riley, 1966)
      real(kind=r8), intent(out), dimension(number_of_records) :: St
      !> total fluoride  (Riley, 1965)
      real(kind=r8), intent(out), dimension(number_of_records) :: Ft
      !> total boron
      !! from either Uppstrom (1974) or Lee et al. (2010), depending on optB
      real(kind=r8), intent(out), dimension(number_of_records) :: Bt

! Local variables
      real(kind=rx) :: ssal
      real(kind=rx) :: p
      real(kind=rx) :: tempot, tempis68, tempot68
      real(kind=rx) :: tempis
      real(kind=r8) :: ionic_strength
      !real(kind=r8) is2,  sqrtis
      real(kind=r8) :: Ks_0p, Kf_0p
      real(kind=r8) :: total2free, free2SWS, total2SWS, SWS2total
      real(kind=r8) :: total2free_0p, free2SWS_0p, total2SWS_0p
! REAL(kind=r8) :: free2SWS, free2SWS_0p

      real(kind=r8) :: dtempot, dtempot68


      !REAL(kind=r8), DIMENSION(12) :: a0, a1, a2, b0, b1, b2
      real(kind=r8), dimension(12) :: lnkpok0

      integer :: i, icount

      real(kind=r8) :: temperature, temperature_kelvin, inverse_temperature_kelvin, log_temperature_kelvin, tk0
      !REAL(kind=r8) :: sqrts, s15, s2
      real(kind=r8) :: salinity, scl

      real(kind=r8) :: atmospheric_pressure, hydrostatic_pressure, total_pressure
! Arrays to pass optional arguments into or use defaults (Dickson et al., 2007)
      character(3) :: opB
      character(2) :: opKf
      character(3) :: opK1K2
      character(7) :: opGAS


! Set defaults for optional arguments (in Fortran 90)
! Note:  Optional arguments with f2py (python) are set above with
!        the !f2py statements that precede each type declaraion
      if (present(optB)) then
         opB = optB
      else
         opB = 'l10'
      endif
      if (present(optKf)) then
         opKf = optKf
      else
         opKf = 'pf'
      endif
      if (present(optK1K2)) then
         opK1K2 = optK1K2
      else
         opK1K2 = 'l'
      endif
      if (present(optGAS)) then
         opGAS = optGAS
      else
         opGAS = 'Pinsitu'
      endif

      icount = 0
      do i = 1, number_of_records
         icount = icount + 1
!    ===============================================================
!    Convert model depth -> press; convert model Theta -> T in situ
!    ===============================================================
!    * Model temperature tracer is usually "potential temperature"
!    * Model vertical grid is usually in meters
!    BUT carbonate chem routines require pressure & in-situ T
!    Thus before computing chemistry, if appropriate,
!    convert these 2 model vars (input to this routine)
!     - depth [m] => convert to pressure [db]
!     - potential temperature (C) => convert to in-situ T (C)
!    -------------------------------------------------------
!    1)  Compute pressure [db] from depth [m] and latitude [degrees] (if input is m, for models)
         if (trim(optP) == 'm' ) then
!       Compute pressure [db] from depth [m] and latitude [degrees]
            p = p80(depth(i), lat(i))
         elseif (trim(optP) == 'db' ) then
!       In this case (where optP = 'db'), p is input & output (no depth->pressure conversion needed)
            p = depth(i)
         else
            print *,"optP must be 'm' or 'db'"
            stop
         endif

!    2) Convert potential T to in-situ T (if input is Tpot, i.e. case for models):
         if (trim(optT) == 'Tpot' .or. trim(optT) == 'tpot') then
            tempot = temp(i)
!       This is the case for most models and some data
!       a) Convert the pot. temp on today's "ITS 90" scale to older IPTS 68 scale
!          (see Dickson et al., Best Practices Guide, 2007, Chap. 5, p. 7, including footnote)
            tempot68 = (tempot - 0.0002_rx) / 0.99975_rx
!       b) Compute "in-situ Temperature" from "Potential Temperature" (both on IPTS 68)
            tempis68 = sw_temp(sal(i), tempot68, p, SGLE(0.0_rx) )
!       c) Convert the in-situ temp on older IPTS 68 scale to modern scale (ITS 90)
            tempis = 0.99975_rx*tempis68 + 0.0002_rx
!       Note: parts (a) and (c) above are tiny corrections;
!             part  (b) is a big correction for deep waters (but zero at surface)
         elseif (trim(optT) == 'Tinsitu' .or. trim(optT) == 'tinsitu') then
!       When optT = 'Tinsitu', tempis is input & output (no tempot needed)
            tempis    = temp(i)
            tempis68  = (temp(i) - 0.0002_rx) / 0.99975_rx
            dtempot68 = sw_ptmp(dble(sal(i)), dble(tempis68), dble(p), 0.0d0)
            dtempot   = 0.99975_rx*dtempot68 + 0.0002_rx
         else
            print *,"optT must be either 'Tpot' or 'Tinsitu'"
            print *,"you specified optT =", trim(optT)
            stop
         endif

!    Compute constants:
         block
            logical :: within_range
            within_range = .false.
            within_range = (temp(i) >= -5.0_rx .and. temp(i) < 1.0e+2_rx)
            if (.not. within_range) then

               k0(i)   = huge(r8)
               k1(i)   = huge(r8)
               k2(i)   = huge(r8)
               Kb(i)   = huge(r8)
               Kw(i)   = huge(r8)
               Ks(i)   = huge(r8)
               Kf(i)   = huge(r8)
               Kspc(i) = huge(r8)
               Kspa(i) = huge(r8)
               K1p(i)  = huge(r8)
               K2p(i)  = huge(r8)
               K3p(i)  = huge(r8)
               Ksi(i)  = huge(r8)
               Bt(i)   = huge(r8)
               Ft(i)   = huge(r8)
               St(i)   = huge(r8)
            else
!       Test to indicate if any of input variables are unreasonable
               if (sal(i) < 0.  .or.  sal(i) > 1e+3) then
                  print *, 'i, icount, temp, sal =', i, icount, temp(i), sal(i)
               endif
!       Zero out negative salinity (prev case for OCMIP2 model w/ slightly negative S in some coastal cells)
               if (sal(i) < 0.0) then
                  ssal = 0.0
               else
                  ssal = sal(i)
               endif

!       Absolute temperature (Kelvin) and related values
               temperature = dble(tempis)
               temperature_kelvin = zero_c_in_kelvin + temperature
               inverse_temperature_kelvin=1.0d0/temperature_kelvin
               log_temperature_kelvin = log(temperature_kelvin)

!       Atmospheric pressure
               atmospheric_pressure = dble(Patm(i))

!       Hydrostatic pressure (prb is in bars)
               hydrostatic_pressure = dble(p) / 10.0d0

!       Salinity and simply related values
               salinity = dble(ssal)

               scl=salinity/knudsen_chlorinity_constant

!       Ionic strength:
               ! more magic numbers
               ionic_strength = 19.924d0*salinity/(1000.0d0 - 1.005d0*salinity)

!       Total concentrations for sulfate, fluoride, and boron

!       Sulfate: Morris & Riley (1966)
               ! TODO, find reference
               St(i) = 0.14d0 * scl/sulfate_atomic_mass

!       Fluoride:  Riley (1965)
               ! atomic mass fluoride
               Ft(i) = 0.000067d0 * scl/fluoride_atomic_mass

!       Boron:
               if (trim(opB) == 'l10') then
!          New formulation from Lee et al (2010)
                  Bt(i) = 0.0002414d0 * scl/boron_atomic_mass
               elseif (trim(opB) == 'u74') then
!          Classic formulation from Uppström (1974)
                  Bt(i) = 0.000232d0  * scl/boron_atomic_mass
               else
                  print *,"optB must be 'l10' or 'u74'"
                  stop
               endif

!       K0 (K Henry)
!       CO2(g) <-> CO2(aq.)
!       K0  = [CO2]/ fCO2
!       Weiss (1974)   [mol/kg/atm]
               total_pressure = calculate_ptot(opgas, atmospheric_pressure, hydrostatic_pressure)

               tk0 = calculate_tk0(opgas, temperature_kelvin, dtempot)

               k0(i) = calculate_k0(tk0, salinity)

!       K1 = [H][HCO3]/[H2CO3]
!       K2 = [H][CO3]/[HCO3]
               k1(i) = calculate_k1(opk1k2, inverse_temperature_kelvin, log_temperature_kelvin, salinity)
               k2(i) = calculate_k2(opk1k2, inverse_temperature_kelvin, log_temperature_kelvin, salinity)

!       Kb = [H][BO2]/[HBO2]
!       (total scale)
!       Millero p.669 (1995) using data from Dickson (1990)
               ! DOI: https://doi.org/10.1016/0016-7037(94)00354-O
               Kb(i) = calculate_kb(temperature_kelvin, inverse_temperature_kelvin, log_temperature_kelvin, salinity)

!       K1p = [H][H2PO4]/[H3PO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
!       Millero (1995), p.670, eq. 65
!       Use Millero equation's 115.540 constant instead of 115.525 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
               k1p(i) = calculate_k1p(inverse_temperature_kelvin, log_temperature_kelvin, salinity)
!       K2p = [H][HPO4]/[H2PO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
!       Millero (1995), p.670, eq. 66
!       Use Millero equation's 172.1033 constant instead of 172.0833 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
               k2p(i) = calculate_k2p(inverse_temperature_kelvin, log_temperature_kelvin, salinity)

!       K3p = [H][PO4]/[HPO4]
!       (seawater scale)
!       DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
!       Millero (1995), p.670, eq. 67
!       Use Millero equation's 18.126 constant instead of 18.141 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
               k3p(i) = calculate_k3p(inverse_temperature_kelvin, salinity)
!       Ksi = [H][SiO(OH)3]/[Si(OH)4]
!       (seawater scale)
!       Millero (1995), p.671, eq. 72
!       Use Millero equation's 117.400 constant instead of 117.385 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
               ksi(i) = calculate_ksi(inverse_temperature_kelvin, log_temperature_kelvin, salinity, ionic_strength)

!       Kw = [H][OH]
!       (seawater scale)
!       Millero (1995) p.670, eq. 63 from composite data
!       Use Millero equation's 148.9802 constant instead of 148.9652 (Dickson et al., 2007).
!       The latter is only an crude approximation to convert to Total scale (by subtracting 0.015)
!       And we want to stay on the SWS scale anyway for the pressure correction later.
               kw(i) = calculate_kw(inverse_temperature_kelvin, log_temperature_kelvin, salinity)


!       Kspc (calcite) - apparent solubility product of calcite
!       (no scale)
!       Kspc = [Ca2+] [CO32-] when soln is in equilibrium w/ calcite
!       Mucci 1983 mol/kg-soln
               kspc(i) = calculate_kspc(temperature_kelvin, salinity)

!       Kspa (aragonite) - apparent solubility product of aragonite
!       (no scale)
!       Kspa = [Ca2+] [CO32-] when soln is in equilibrium w/ aragonite
!       Mucci 1983 mol/kg-soln
               kspa(i) = calculate_kspa(temperature_kelvin, salinity)

!       Ks = [H][SO4]/[HSO4]
!       (free scale)
!       Dickson (1990, J. chem. Thermodynamics 22, 113)
               Ks_0p = calculate_ks_no_pressure(inverse_temperature_kelvin, log_temperature_kelvin, salinity, ionic_strength)

!       Kf = [H][F]/[HF]
!       (total scale)
               kf_0p = calculate_kf_no_pressure(opkf, inverse_temperature_kelvin, ionic_strength, salinity, st(i), ks_0p)
!       Pressure effect on all other K's (based on Millero, (1995)
!           index: K1(1), K2(2), Kb(3), Kw(4), Ks(5), Kf(6), Kspc(7), Kspa(8),
!                  K1p(9), K2p(10), K3p(11), Ksi(12)
               call calculate_all_pressure_correction_factors(temperature, hydrostatic_pressure, temperature_kelvin, lnkpok0)

!       Pressure effect on K0 based on Weiss (1974, equation 5)
               k0(i) = k0(i) * exp( ((1-total_pressure)*co2_partial_molar_volume)/(ideal_gas_constant_codata*tk0) )

!       Pressure correction on Ks (Free scale)
               Ks(i) = Ks_0p*exp(lnkpok0(5))
!       Conversion factor total -> free scale
               total2free     = 1.d0/(1.d0 + St(i)/Ks(i))   ! Kfree = Ktotal*total2free
!       Conversion factor total -> free scale at pressure zero
               total2free_0p  = 1.d0/(1.d0 + St(i)/Ks_0p)   ! Kfree = Ktotal*total2free

!       Pressure correction on Kf
!       Kf must be on FREE scale before correction
               Kf_0p = Kf_0p * total2free_0p   !Convert from Total to Free scale (pressure 0)
               Kf(i) = Kf_0p * exp(lnkpok0(6)) !Pressure correction (on Free scale)
               Kf(i) = Kf(i)/total2free        !Convert back from Free to Total scale

!       Convert between seawater and total hydrogen (pH) scales
               free2SWS  = 1.d0 + St(i)/Ks(i) + Ft(i)/(Kf(i)*total2free)  ! using Kf on free scale
               total2SWS = total2free * free2SWS                          ! KSWS = Ktotal*total2SWS
               SWS2total = 1.d0 / total2SWS

!       Conversion at pressure zero
               free2SWS_0p  = 1.d0 + St(i)/Ks_0p + Ft(i)/(Kf_0p)  ! using Kf on free scale
               total2SWS_0p = total2free_0p * free2SWS_0p         ! KSWS = Ktotal*total2SWS

!       Convert from Total to Seawater scale before pressure correction
!       Must change to SEAWATER scale: K1, K2, Kb
               if (trim(optK1K2) == 'l') then
                  k1(i)  = k1(i)*total2SWS_0p
                  k2(i)  = k2(i)*total2SWS_0p
                  !This conversion is unnecessary for the K1,K2 from Millero (2010),
                  !since we use here the formulation already on the seawater scale
               endif
               Kb(i)  = Kb(i)*total2SWS_0p

!       Already on SEAWATER scale: K1p, K2p, K3p, Kb, Ksi, Kw

!       Other contants (keep on another scale):
!          - K0         (independent of pH scale, already pressure corrected)
!          - Ks         (already on Free scale;   already pressure corrected)
!          - Kf         (already on Total scale;  already pressure corrected)
!          - Kspc, Kspa (independent of pH scale; pressure-corrected below)

!       Perform actual pressure correction (on seawater scale)
               k1(i)   = k1(i)*exp(lnkpok0(1))
               k2(i)   = k2(i)*exp(lnkpok0(2))
               Kb(i)   = Kb(i)*exp(lnkpok0(3))
               Kw(i)   = Kw(i)*exp(lnkpok0(4))
               Kspc(i) = Kspc(i)*exp(lnkpok0(7))
               Kspa(i) = Kspa(i)*exp(lnkpok0(8))
               K1p(i)  = K1p(i)*exp(lnkpok0(9))
               K2p(i)  = K2p(i)*exp(lnkpok0(10))
               K3p(i)  = K3p(i)*exp(lnkpok0(11))
               Ksi(i)  = Ksi(i)*exp(lnkpok0(12))

!       Convert back to original total scale:
               k1(i)  = k1(i) *SWS2total
               k2(i)  = k2(i) *SWS2total
               K1p(i) = K1p(i)*SWS2total
               K2p(i) = K2p(i)*SWS2total
               K3p(i) = K3p(i)*SWS2total
               Kb(i)  = Kb(i) *SWS2total
               Ksi(i) = Ksi(i)*SWS2total
               Kw(i)  = Kw(i) *SWS2total




            endif
         end block

      end do

      return
   end subroutine constants

   function calculate_tk0(op_gas, tk, dtempot) result(tk0_value)
      ! Arguments
      character(len=*), intent(in) :: op_gas
      real(r8), intent(in) :: tk          ! in situ temperature (K)
      real(r8), intent(in) :: dtempot     ! potential temperature (C)

      ! Result
      real(r8) :: tk0_value

      ! Determine temperature based on gas option
      select case (trim(adjustl(op_gas)))
       case ('Pzero', 'pzero')
         tk0_value = tk                        ! in situ temperature (K)
       case ('Ppot', 'ppot')
         tk0_value = dtempot + zero_c_in_kelvin ! potential temperature (K)
       case ('Pinsitu', 'pinsitu')
         tk0_value = tk                         ! in situ temperature (K)
       case default
         print *, "error: op_gas must be 'Pzero'/'pzero', 'Ppot'/'ppot', or 'Pinsitu'/'pinsitu'"
         stop
      end select


   end function calculate_tk0

   function calculate_ptot(op_gas, patmd, prb) result(ptot)
      ! Arguments
      character(len=*), intent(in) :: op_gas
      real(r8), intent(in) :: patmd       ! atmospheric pressure (atm)
      real(r8), intent(in) :: prb         ! hydrostatic pressure (bar)

      ! Result
      real(r8) :: ptot

      ! Local variables
      real(r8) :: phydro_atm

      ! Determine pressure based on gas option
      select case (trim(adjustl(op_gas)))
       case ('Pzero', 'pzero')
         ptot = patmd                    ! atmospheric pressure only

       case ('Ppot', 'ppot')
         ptot = patmd                     ! atmospheric pressure only

       case ('Pinsitu', 'pinsitu')
         phydro_atm = prb / bar_to_atm    ! convert bar to atm
         ptot = patmd + phydro_atm        ! total pressure

       case default
         print *, "error: op_gas must be 'Pzero'/'pzero', 'Ppot'/'ppot', or 'Pinsitu'/'pinsitu'"
         stop
      end select


   end function calculate_ptot

   function calculate_k0(tk0, s) result(k0_value)
      ! Arguments
      real(r8), intent(in) :: tk0          ! in situ temperature (K)
      real(r8), intent(in) :: s           ! salinity

      ! Result
      real(r8) :: k0_value

      ! Local variables
      real(r8) :: tmp, nk0we74

      ! Calculate K0 (Weiss 1974 formulation)
      tmp = 9345.17_r8/tk0 - 60.2409_r8 + 23.3585_r8 * log(tk0/100.0_r8)
      nk0we74 = tmp + s*(0.023517_r8 - 0.00023656_r8*tk0 + 0.0047036e-4_r8*tk0*tk0)
      k0_value = exp(nk0we74)

   end function calculate_k0

   function calculate_k1(op_k1k2, invtk, dlogtk, s) result(k1_value)

      ! Arguments
      character(len=*), intent(in) :: op_k1k2
      real(r8), intent(in) :: invtk, dlogtk, s
      real(r8) :: s2, sqrts

      ! Result
      real(r8) :: k1_value

      ! Local variables
      real(r8) :: pk1o, ma1, mb1, mc1, pk1

      s2 = s*s
      sqrts = sqrt(s)

      select case (trim(op_k1k2))
       case ('l')
         ! Mehrbach et al. (1973) refit, by Lueker et al. (2000) (total scale)
         k1_value = 10.0_r8**(-1.0_r8*(3633.86_r8*invtk - 61.2172_r8 + 9.6777_r8*dlogtk &
            - 0.011555_r8*s + 0.0001152_r8*s2))

       case ('m10')
         ! Millero (2010, Mar. Fresh Wat. Res.) (seawater scale)
         pk1o = 6320.813_r8*invtk + 19.568224_r8*dlogtk - 126.34048_r8
         ma1 = 13.4038_r8*sqrts + 0.03206_r8*s - (5.242e-5_r8)*s2
         mb1 = -530.659_r8*sqrts - 5.8210_r8*s
         mc1 = -2.0664_r8*sqrts
         pk1 = pk1o + ma1 + mb1*invtk + mc1*dlogtk
         k1_value = 10.0_r8**(-pk1)

       case ('w14')
         ! Waters, Millero, Woosley (Mar. Chem., 165, 66-67, 2014) (seawater scale)
         pk1o = 6320.813_r8*invtk + 19.568224_r8*dlogtk - 126.34048_r8
         ma1 = 13.409160_r8*sqrts + 0.031646_r8*s - (5.1895e-5_r8)*s2
         mb1 = -531.3642_r8*sqrts - 5.713_r8*s
         mc1 = -2.0669166_r8*sqrts
         pk1 = pk1o + ma1 + mb1*invtk + mc1*dlogtk
         k1_value = 10.0_r8**(-pk1)

       case default
         print *, "error: op_k1k2 must be 'l', 'm10', or 'w14'"
         stop
      end select

   end function calculate_k1

   function calculate_k2(op_k1k2, invtk, dlogtk, s) result(k2_value)
      implicit none

      ! Arguments
      character(len=*), intent(in) :: op_k1k2
      real(r8), intent(in) :: invtk, dlogtk, s
      real(r8) :: s2, sqrts


      ! Result
      real(r8) :: k2_value

      ! Local variables
      real(r8) :: pk2o, ma2, mb2, mc2, pk2

      s2 = s*s
      sqrts = sqrt(s)

      select case (trim(op_k1k2))
       case ('l')
         ! Mehrbach et al. (1973) refit, by Lueker et al. (2000) (total scale)
         k2_value = 10.0_r8**(-1.0_r8*(471.78_r8*invtk + 25.9290_r8 - 3.16967_r8*dlogtk &
            - 0.01781_r8*s + 0.0001122_r8*s2))

       case ('m10')
         ! Millero (2010, Mar. Fresh Wat. Res.) (seawater scale)
         pk2o = 5143.692_r8*invtk + 14.613358_r8*dlogtk - 90.18333_r8
         ma2 = 21.3728_r8*sqrts + 0.1218_r8*s - (3.688e-4_r8)*s2
         mb2 = -788.289_r8*sqrts - 19.189_r8*s
         mc2 = -3.374_r8*sqrts
         pk2 = pk2o + ma2 + mb2*invtk + mc2*dlogtk
         k2_value = 10.0_r8**(-pk2)

       case ('w14')
         ! Waters, Millero, Woosley (Mar. Chem., 165, 66-67, 2014) (seawater scale)
         pk2o = 5143.692_r8*invtk + 14.613358_r8*dlogtk - 90.18333_r8
         ma2 = 21.225890_r8*sqrts + 0.12450870_r8*s - (3.7243e-4_r8)*s2
         mb2 = -779.3444_r8*sqrts - 19.91739_r8*s
         mc2 = -3.3534679_r8*sqrts
         pk2 = pk2o + ma2 + mb2*invtk + mc2*dlogtk
         k2_value = 10.0_r8**(-pk2)

       case default
         print *, "error: op_k1k2 must be 'l', 'm10', or 'w14'"
         stop
      end select

   end function calculate_k2

   function calculate_kb(tk, invtk, dlogtk, s) result(kb_value)

      ! Arguments
      real(r8), intent(in) :: tk, invtk, dlogtk, s
      real(r8) :: sqrts, s15, s2

      ! Result
      real(r8) :: kb_value

      sqrts = sqrt(s)
      s15 = s**1.5_r8
      s2 = s*s
      ! Calculate Kb (total scale)
      kb_value = exp((-8966.90_r8 - 2890.53_r8*sqrts - 77.942_r8*s + &
         1.728_r8*s15 - 0.0996_r8*s2)*invtk + &
         (148.0248_r8 + 137.1942_r8*sqrts + 1.62142_r8*s) + &
         (-24.4344_r8 - 25.085_r8*sqrts - 0.2474_r8*s) * &
         dlogtk + 0.053105_r8*sqrts*tk)


   end function calculate_kb

   function calculate_k1p(invtk, dlogtk, s) result(k1p_value)

      ! Arguments
      real(r8), intent(in) :: invtk, dlogtk, s
      real(r8) :: sqrts

      ! Result
      real(r8) :: k1p_value

      sqrts = sqrt(s)

      ! K1p = [H][H2PO4]/[H3PO4] (seawater scale)
      ! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
      k1p_value = exp(-4576.752_r8*invtk + 115.540_r8 - 18.453_r8*dlogtk + &
         (-106.736_r8*invtk + 0.69171_r8) * sqrts + &
         (-0.65643_r8*invtk - 0.01844_r8) * s)

   end function calculate_k1p

   function calculate_k2p(invtk, dlogtk, s) result(k2p_value)

      ! Arguments
      real(r8), intent(in) :: invtk, dlogtk, s
      real(r8) :: sqrts

      ! Result
      real(r8) :: k2p_value

      sqrts = sqrt(s)
      ! K2p = [H][HPO4]/[H2PO4] (seawater scale)
      ! DOE(1994) eq 7.2.23 with footnote using data from Millero (1974)
      k2p_value = exp(-8814.715_r8*invtk + 172.1033_r8 - 27.927_r8*dlogtk + &
         (-160.340_r8*invtk + 1.3566_r8)*sqrts + &
         (0.37335_r8*invtk - 0.05778_r8)*s)

   end function calculate_k2p

   function calculate_k3p(invtk, s) result(k3p_value)

      ! Arguments
      real(r8), intent(in) :: invtk, s
      real(r8) :: sqrts

      ! Result
      real(r8) :: k3p_value
      sqrts = sqrt(s)
      ! K3p = [H][PO4]/[HPO4] (seawater scale)
      ! DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
      k3p_value = exp(-3070.75_r8*invtk - 18.126_r8 + &
         (17.27039_r8*invtk + 2.81197_r8) * sqrts + &
         (-44.99486_r8*invtk - 0.09984_r8) * s)

   end function calculate_k3p

   function calculate_ksi(invtk, dlogtk, s, is) result(ksi_value)

      ! Arguments
      real(r8), intent(in) :: invtk, dlogtk, is, s
      real(r8) :: sqrtis, is2

      ! Result
      real(r8) :: ksi_value

      is2 = is * is
      sqrtis = sqrt(is)

      ! Ksi = [H][SiO(OH)3]/[Si(OH)4] (seawater scale)
      ! Millero (1995), p.671, eq. 72
      ksi_value = exp(-8904.2_r8*invtk + 117.400_r8 - 19.334_r8*dlogtk + &
         (-458.79_r8*invtk + 3.5913_r8) * sqrtis + &
         (188.74_r8*invtk - 1.5998_r8) * is + &
         (-12.1652_r8*invtk + 0.07871_r8) * is2 + &
         log(1.0_r8 - 0.001005_r8*s))

   end function calculate_ksi

   function calculate_kw(invtk, dlogtk, s) result(kw_value)

      ! Arguments
      real(r8), intent(in) :: invtk, dlogtk, s
      real(r8) :: sqrts

      ! Result
      real(r8) :: kw_value

      sqrts = sqrt(s)

      ! Kw = [H][OH] (seawater scale)
      ! Millero (1995) p.670, eq. 63 from composite data
      kw_value = exp(-13847.26_r8*invtk + 148.9802_r8 - 23.6521_r8*dlogtk + &
         (118.67_r8*invtk - 5.977_r8 + 1.0495_r8 * dlogtk) * sqrts - &
         0.01615_r8 * s)

   end function calculate_kw

   function calculate_kspc(tk, s) result(kspc_value)

      ! Arguments
      real(r8), intent(in) :: tk, s
      real(r8) :: sqrts, s15

      ! Result
      real(r8) :: kspc_value

      s15 = s**1.5_r8
      sqrts = sqrt(s)
      ! Kspc (calcite) - apparent solubility product of calcite
      ! Kspc = [Ca2+] [CO32-] when soln is in equilibrium w/ calcite
      ! Mucci 1983 mol/kg-soln
      kspc_value = 10.0_r8**(-171.9065_r8 - 0.077993_r8*tk + 2839.319_r8/tk + &
         71.595_r8*log10(tk) + &
         (-0.77712_r8 + 0.0028426_r8*tk + 178.34_r8/tk)*sqrts - &
         0.07711_r8*s + 0.0041249_r8*s15)

   end function calculate_kspc

   function calculate_kspa(tk, s) result(kspa_value)

      ! Arguments
      real(r8), intent(in) :: tk, s
      real(r8) :: sqrts, s15

      ! Result
      real(r8) :: kspa_value
      s15 = s**1.5_r8
      sqrts = sqrt(s)

      ! Kspa (aragonite) - apparent solubility product of aragonite
      ! Kspa = [Ca2+] [CO32-] when soln is in equilibrium w/ aragonite
      ! Mucci 1983 mol/kg-soln
      kspa_value = 10.0_r8**(-171.945_r8 - 0.077993_r8*tk + 2903.293_r8/tk + &
         71.595_r8*log10(tk) + &
         (-0.068393_r8 + 0.0017276_r8*tk + 88.135_r8/tk)*sqrts - &
         0.10018_r8*s + 0.0059415_r8*s15)

   end function calculate_kspa

   function calculate_ks_no_pressure(invtk, dlogtk, s, is) result(ks_0p_value)

      ! Arguments
      real(r8), intent(in) :: invtk, dlogtk, s, is
      real(r8) :: sqrtis, is2

      ! Result
      real(r8) :: ks_0p_value
      is2 = is * is
      sqrtis = sqrt(is)

      ! Ks = [H][SO4]/[HSO4] (free scale) at zero pressure
      ! Dickson (1990, J. chem. Thermodynamics 22, 113)
      ks_0p_value =  exp(-4276.1_r8*invtk + 141.328_r8 - 23.093_r8*dlogtk + &
         (-13856._r8*invtk + 324.57_r8 - 47.986_r8*dlogtk) * sqrtis + &
         (35474._r8*invtk - 771.54_r8 + 114.723_r8*dlogtk) * is &
         - 2698._r8*invtk*is**1.5_r8  + 1776._r8*invtk*is2 + &
         log(1.0_r8 - 0.001005_r8*s))

   end function calculate_ks_no_pressure

   function calculate_kf_no_pressure(op_kf, invtk, is, s, st, ks_0p) result(kf_0p_value)
      implicit none

      ! Arguments
      character(len=*), intent(in) :: op_kf
      real(r8), intent(in) :: invtk, is, s, st, ks_0p
      real(r8) :: sqrtis, sqrts


      ! Result
      real(r8) :: kf_0p_value

      sqrtis = sqrt(is)
      sqrts = sqrt(s)
      ! Kf = [H][F]/[HF] (total scale) at zero pressure
      select case (trim(op_kf))
       case ('dg')
         ! Dickson and Riley (1979) -- change pH scale to total (following Dickson & Goyet, 1994)
         kf_0p_value = exp(1590.2_r8*invtk - 12.641_r8 + 1.525_r8*sqrtis + &
            log(1.0_r8 - 0.001005_r8*s) + &
            log(1.0_r8 + st/ks_0p))

       case ('pf')
         ! Perez and Fraga (1987) - Already on Total scale
         ! Formulation as given in Dickson et al. (2007)
         kf_0p_value = exp(874._r8*invtk - 9.68_r8 + 0.111_r8*sqrts)

       case default
         print *, "error: op_kf must be either 'dg' or 'pf'"
         stop
      end select

   end function calculate_kf_no_pressure

   subroutine calculate_all_pressure_correction_factors(t, prb, tk, pressure_correction_factors)

      ! Arguments
      real(r8), intent(in) :: t, prb, tk
      real(r8), intent(out) :: pressure_correction_factors(12)


      ! Local variables
      integer :: ipc
      real(r8) :: deltav, deltak
      do ipc = 1, 12
         deltav = a0(ipc) + a1(ipc)*t + a2(ipc)*t*t
         deltak = b0(ipc) + b1(ipc)*t + b2(ipc)*t*t
         pressure_correction_factors(ipc) = (-(deltav) + (0.5_r8*deltak*prb)) * prb / (R_jkmol_scaled_by_10*tk)
      end do

   end subroutine calculate_all_pressure_correction_factors



end module mocsy_constants
