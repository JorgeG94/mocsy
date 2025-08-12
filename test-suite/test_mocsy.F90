module test_mocsy
  use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
  use, intrinsic :: iso_fortran_env, only: error_unit
   USE mocsy_singledouble
   USE mocsy_constants
   USE mocsy_vars
   USE mocsy_derivauto
   use mocsy_test_helpers, only: is_equal
  implicit none 

  public :: collect_mocsy_suite


contains 

subroutine collect_mocsy_suite(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ & 
    new_unittest("constants", test_constants), & 
    new_unittest("kzero", test_kzero), &
    new_unittest("derivauto", test_derivauto), &
    new_unittest("derivnum", test_derivnum), &
    new_unittest("errors", test_errors), &
    new_unittest("kprime", test_kprime), &
    new_unittest("phizero", test_phizero), &
    new_unittest("buffesm", test_buffesm) &
  ]

end subroutine collect_mocsy_suite

subroutine test_constants(error) 
type(error_type), allocatable, intent(out) :: error
!  For vars routine (called below)
!  "vars" Output variables:
   REAL(kind=rx), DIMENSION(100) :: ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis
!  "vars" Input variables
   INTEGER :: N
   REAL(kind=rx), DIMENSION(100) :: temp, sal, alk, dic, sil, phos, Patm, depth, lat
   REAL(kind=rx), DIMENSION(6,100) :: ph_deriv, pco2_deriv, fco2_deriv, co2_deriv, hco3_deriv, co3_deriv, OmegaA_deriv, OmegaC_deriv
   REAL(kind=rx) ::  gamma_DIC, gamma_Alk, beta_DIC, beta_Alk, omega_DIC, omega_Alk
!  "vars" Input options
   CHARACTER(10) :: optCON, optT, optP, optB, optKf, optK1K2,optGAS

  !> solubility of CO2 in seawater (Weiss, 1974), also known as K0
  REAL(kind=r8), DIMENSION(6) :: K0
  !> K1 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=r8), DIMENSION(6) :: K1
  !> K2 for the dissociation of carbonic acid from Lueker et al. (2000) or Millero (2010), depending on optK1K2
  REAL(kind=r8), DIMENSION(6) :: K2
  !> equilibrium constant for dissociation of boric acid 
  REAL(kind=r8), DIMENSION(6) :: Kb
  !> equilibrium constant for the dissociation of water (Millero, 1995)
  REAL(kind=r8), DIMENSION(6) :: Kw
  !> equilibrium constant for the dissociation of bisulfate (Dickson, 1990)
  REAL(kind=r8), DIMENSION(6) :: Ks
  !> equilibrium constant for the dissociation of hydrogen fluoride 
  !! either from Dickson and Riley (1979) or from Perez and Fraga (1987), depending on optKf
  REAL(kind=r8), DIMENSION(6) :: Kf
  !> solubility product for calcite (Mucci, 1983)
  REAL(kind=r8), DIMENSION(6) :: Kspc
  !> solubility product for aragonite (Mucci, 1983)
  REAL(kind=r8), DIMENSION(6) :: Kspa
  !> 1st dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), DIMENSION(6) :: K1p
  !> 2nd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), DIMENSION(6) :: K2p
  !> 3rd dissociation constant for phosphoric acid (Millero, 1995)
  REAL(kind=r8), DIMENSION(6) :: K3p
  !> equilibrium constant for the dissociation of silicic acid (Millero, 1995)
  REAL(kind=r8), DIMENSION(6) :: Ksi
  !> total sulfate (Morris & Riley, 1966)
  REAL(kind=r8), DIMENSION(6) :: St
  !> total fluoride  (Riley, 1965)
  REAL(kind=r8), DIMENSION(6) :: Ft
  !> total boron
  !! from either Uppstrom (1974) or Lee et al. (2010), depending on optB
  REAL(kind=r8), DIMENSION(6) :: Bt

! Parameter arrays for correctness checking
real(kind=r8), parameter :: K0_ref(6) = [ &
    0.058223497769_r8, 0.050556962375_r8, 0.043899912278_r8, &
    0.038119424259_r8, 0.033100077664_r8, 0.028741649767_r8 ]

real(kind=r8), parameter :: K1_ref(6) = [ &
    0.000000814647_r8, 0.000000909626_r8, 0.000001014354_r8, &
    0.000001129668_r8, 0.000001256455_r8, 0.000001395656_r8 ]

real(kind=r8), parameter :: K2_ref(6) = [ &
    0.000000000444_r8, 0.000000000476_r8, 0.000000000511_r8, &
    0.000000000548_r8, 0.000000000589_r8, 0.000000000632_r8 ]

real(kind=r8), parameter :: Kb_ref(6) = [ &
    0.000000001303_r8, 0.000000001480_r8, 0.000000001679_r8, &
    0.000000001902_r8, 0.000000002152_r8, 0.000000002432_r8 ]

real(kind=r8), parameter :: Kw_ref(6) = [ &
    0.000000000000_r8, 0.000000000000_r8, 0.000000000000_r8, &
    0.000000000000_r8, 0.000000000000_r8, 0.000000000000_r8 ]

real(kind=r8), parameter :: Ks_ref(6) = [ &
    0.260528321264_r8, 0.281507695490_r8, 0.303598633484_r8, &
    0.326801140739_r8, 0.351108647315_r8, 0.376507557419_r8 ]

real(kind=r8), parameter :: Kf_ref(6) = [ &
    0.004094821803_r8, 0.004239407350_r8, 0.004384697209_r8, &
    0.004530132806_r8, 0.004675137446_r8, 0.004819119223_r8 ]

real(kind=r8), parameter :: Kspc_ref(6) = [ &
    0.000000429916_r8, 0.000000528309_r8, 0.000000646101_r8, &
    0.000000786357_r8, 0.000000952462_r8, 0.000001148108_r8 ]

real(kind=r8), parameter :: Kspa_ref(6) = [ &
    0.000000683037_r8, 0.000000829150_r8, 0.000001001682_r8, &
    0.000001204298_r8, 0.000001440941_r8, 0.000001715797_r8 ]

real(kind=r8), parameter :: K1p_ref(6) = [ &
    0.024768086450_r8, 0.026361992774_r8, 0.028025832471_r8, &
    0.029760061616_r8, 0.031564919496_r8, 0.033440413999_r8 ]

real(kind=r8), parameter :: K2p_ref(6) = [ &
    0.000000666884_r8, 0.000000736317_r8, 0.000000811186_r8, &
    0.000000891700_r8, 0.000000978047_r8, 0.000001070394_r8 ]

real(kind=r8), parameter :: K3p_ref(6) = [ &
    0.000000000454_r8, 0.000000000509_r8, 0.000000000569_r8, &
    0.000000000636_r8, 0.000000000709_r8, 0.000000000789_r8 ]

real(kind=r8), parameter :: Ksi_ref(6) = [ &
    0.000000000149_r8, 0.000000000169_r8, 0.000000000192_r8, &
    0.000000000218_r8, 0.000000000247_r8, 0.000000000279_r8 ]

real(kind=r8), parameter :: St_ref(6) = [ &
    0.028235434133_r8, 0.028235434133_r8, 0.028235434133_r8, &
    0.028235434133_r8, 0.028235434133_r8, 0.028235434133_r8 ]

real(kind=r8), parameter :: Ft_ref(6) = [ &
    0.000068324401_r8, 0.000068324401_r8, 0.000068324401_r8, &
    0.000068324401_r8, 0.000068324401_r8, 0.000068324401_r8 ]

real(kind=r8), parameter :: Bt_ref(6) = [ &
    0.000432602930_r8, 0.000432602930_r8, 0.000432602930_r8, &
    0.000432602930_r8, 0.000432602930_r8, 0.000432602930_r8 ]


!  Local variables:
   INTEGER :: i


!> Typical options for observations
   optCON  = 'mol/kg'  ! input concentrations are in MOL/KG
   optT    = 'Tinsitu' ! input temperature, variable 'temp' is actually IN SITU temp [°C]
   optP    = 'db'      ! input variable 'depth' is in 'DECIBARS'
   optB    = 'l10'
   optK1K2 = 'l'
   optKf   = 'dg'
   optGAS  = 'Pzero'

   DO i = 1,6
     temp(i)   = 2.0            !Can be "Potential temperature" or "In situ temperature" (see optT below)
     sal(i)    = 35.0           !Salinity (practical scale)
     alk(i)    = 2295.*1.e-6      ! Convert obs. S. Ocean ave surf ALK (umol/kg) to mocsy data units (mol/kg)
     dic(i)    = 2154.*1.e-6      ! Convert obs. S. Ocean ave surf DIC (umol/kg) to mocsy data units (mol/kg)
     sil(i)    = 0.
     phos(i)   = 0.
     depth(i) = real(i-1) * 1000. ! Vary depth from 0 to 5000 db by 1000 db
     Patm(i)   = 1.0            !Atmospheric pressure (atm)
     N = i
   END DO


!  Need to change to call latest version of constants & derivnum (to get constand
  call constants (K0, K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa,                   &
                     K1p, K2p, K3p, Ksi,                                    &
                     St, Ft, Bt,                                            &
                     temp, sal, Patm,                                       &
                     depth, lat, 6,                                         &
                     optT='Tinsitu', optP='db', optB='l10', optK1K2=optK1K2, optKf='dg',  &
                     optGAS='Pinsitu')

  call check(error, all(is_equal(K0, K0_ref)), .true., "K0 does not match!")
  if(allocated(error)) return

  call check(error, all(is_equal(K1, K1_ref)), .true., "K1 does not match!")
  if(allocated(error)) return

  call check(error, all(is_equal(K2, K2_ref)), .true., "K2 does not match!")
  if(allocated(error)) return

  call check(error, all(is_equal(Kb, Kb_ref)), .true., "Kb does not match!")
  if(allocated(error)) return

  call check(error, all(is_equal(Kw, Kw_ref)), .true., "Kw does not match!")
  if(allocated(error)) return

  call check(error, all(is_equal(Ks, Ks_ref)), .true., "Ks does not match!")
  if(allocated(error)) return

  call check(error, all(is_equal(Kf, Kf_ref)), .true., "Kf does not match!")
  if(allocated(error)) return

  call check(error, all(is_equal(Kspc, Kspc_ref)), .true., "Kspc does not match!")
  if(allocated(error)) return

  call check(error, all(is_equal(Kspa, Kspa_ref)), .true., "Kspa does not match!")
  if(allocated(error)) return

  call check(error, all(is_equal(K1p, K1p_ref)), .true., "K1p does not match!")
  if(allocated(error)) return

  call check(error, all(is_equal(K2p, K2p_ref)), .true., "K2p does not match!")
  if(allocated(error)) return

  call check(error, all(is_equal(K3p, K3p_ref)), .true., "K3p does not match!")
  if(allocated(error)) return

  call check(error, all(is_equal(Ksi, Ksi_ref)), .true., "Ksi does not match!")
  if(allocated(error)) return

  call check(error, all(is_equal(St, St_ref)), .true., "St does not match!")
  if(allocated(error)) return

  call check(error, all(is_equal(Ft, Ft_ref)), .true., "Ft does not match!")
  if(allocated(error)) return

  call check(error, all(is_equal(Bt, Bt_ref)), .true., "Bt does not match!")
  if(allocated(error)) return


end subroutine test_constants 


subroutine test_kzero(error)
  use mocsy_gasx
type(error_type), allocatable, intent(out) :: error
   INTEGER, PARAMETER :: n = 1
   real(r8), parameter :: k0_co2_ref   =    5.5850871408260290E-002_r8
   real(r8), parameter :: k0_n2o_ref   =    4.1320896614892258E-002_r8

!  Computed variables:
   REAL(kind=r8), DIMENSION(1) :: k0_co2, k0_n2o

!  Input variables
   REAL(kind=rx), DIMENSION(1) :: temp, sal

!  Input at standard T and S
!  temp(1)   = 25.0
!  sal(1)    = 35.0

   temp(1)   = 4.0
   sal(1)    = 34.0

   call kzero('co2',   temp, sal, n, k0_co2)
   call kzero('n2o',   temp, sal, n, k0_n2o)

   call check(error, is_equal(k0_co2_ref, k0_co2(1)), .true., "K0 is not equal!")
   if(allocated(error)) return
   call check(error, is_equal(k0_n2o_ref, k0_n2o(1)), .true., "K0 is not equal!")
   if(allocated(error)) return



end subroutine test_kzero

subroutine test_derivauto(error)
    type(error_type), allocatable, intent(out) :: error
    ! Parameter arrays for derivative correctness checking
! Parameter arrays for derivative correctness checking
real(r8), parameter :: H_deriv_ref(6) = [ &
    -0.000026065467_r8, 0.000028597336_r8, 0.000029356213_r8, &
    0.000001069695_r8, 0.000000000257_r8, 0.000000000211_r8 ]

real(r8), parameter :: pco2_deriv_ref(6) = [ &
    -1214835.394334753510_r8, 1487033.383599504596_r8, 1368207.497009294108_r8, &
    49855.356786936441_r8, 12.826699956291_r8, 8.356157606934_r8 ]

real(r8), parameter :: fco2_deriv_ref(6) = [ &
    -1210608.266623582691_r8, 1481859.118795759976_r8, 1363446.697437425610_r8, &
    49681.880642590346_r8, 12.795489835497_r8, 8.327081613969_r8 ]

real(r8), parameter :: co2_deriv_ref(6) = [ &
    -0.041508862061_r8, 0.050809405033_r8, 0.046749326311_r8, &
    0.001703472863_r8, 0.000000136219_r8, 0.000000227711_r8 ]

real(r8), parameter :: hco3_deriv_ref(6) = [ &
    -0.622426394027_r8, 1.575213620727_r8, 0.701007282646_r8, &
    0.025543616932_r8, -0.000000574614_r8, 0.000000934866_r8 ]

real(r8), parameter :: co3_deriv_ref(6) = [ &
    0.663934263487_r8, -0.626021936743_r8, -0.747755491041_r8, &
    -0.027247049060_r8, 0.000000438429_r8, -0.000001162563_r8 ]

real(r8), parameter :: omegaa_deriv_ref(6) = [ &
    10253.261022613788_r8, -9667.773869045817_r8, -11547.727918822662_r8, &
    -420.781275303493_r8, 0.016808891616_r8, -0.035806210176_r8 ]

real(r8), parameter :: omegac_deriv_ref(6) = [ &
    15855.789129259514_r8, -14950.383441801694_r8, -17857.571205793713_r8, &
    -650.702167440996_r8, 0.014105012189_r8, -0.062941416251_r8 ]


!  For vars routine (called below)
!  "vars" Output variables:
   REAL(kind=rx), DIMENSION(1) :: ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis
!  "vars" Input variables
   INTEGER :: N
   REAL(kind=rx), DIMENSION(1) :: temp, sal, alk, dic, sil, phos, Patm, depth, lat
   REAL(kind=rx), DIMENSION(6,1) :: ph_deriv, pco2_deriv, fco2_deriv, co2_deriv, &
                     hco3_deriv, co3_deriv, omegaa_deriv, omegac_deriv
!  "vars" Input options
   CHARACTER(10) :: optCON, optT, optP, optB, optKf, optK1K2


!  Local variables:
   INTEGER :: i
   REAL(kind=r8) :: H
   REAL(kind=r8), DIMENSION(6) :: H_deriv
   CHARACTER*4 :: invar(6)
   
   invar(1) = 'Alk '
   invar(2) = 'DIC '
   invar(3) = 'Phos'
   invar(4) = 'Sil '
   invar(5) = 'T   '
   invar(6) = 'S   '

!> Typical options for observations
   optCON  = 'mol/kg'  ! input concentrations are in MOL/KG
   optT    = 'Tinsitu' ! input temperature, variable 'temp' is actually IN SITU temp [°C]
   optP    = 'm'      ! input variable 'depth' is in 'DECIBARS'
   optB    = 'l10'
   optK1K2 = 'l'
   optKf   = 'dg'
!> Simple input data (with CONCENTRATION units typical for DATA)
!> (based on observed average surface concentrations from S. Ocean (south of 60°S)--GLODAP and WOA2009)
   DO i = 1,1
     temp(i)   = 18.0            !Can be "Potential temperature" or "In situ temperature" (see optT below)
     sal(i)    = 35.0           !Salinity (practical scale)
     alk(i)    = 2300.*1.e-6      ! Convert obs. S. Ocean ave surf ALK (umol/kg) to mocsy data units (mol/kg)
     dic(i)    = 2000.*1.e-6      ! Convert obs. S. Ocean ave surf DIC (umol/kg) to mocsy data units (mol/kg)
     sil(i)    = 60.*1.e-6
     phos(i)   = 2.*1.e-6
     depth(i)  = 0.
     Patm(i)   = 1.0            !Atmospheric pressure (atm)
     lat(i)    = 0.
     N = i
   END DO

   call vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,         &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,                             &  ! INPUT
             optCON, optT, optP, optB=optB, optK1K2=optK1K2, optKf=optKf    )

   call derivauto(ph_deriv, pco2_deriv, fco2_deriv, co2_deriv, hco3_deriv, co3_deriv,   &
                omegaa_deriv, omegac_deriv,                                             &
                temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,                    &
                optCON, optT, optP, optB=optB, optK1K2=optK1K2, optKf=optKf )

   !call set_precision(20)

    ! [H+] concentration
    H = 10.0**(-ph(1))
    ! derivative of [H+] deduced from that of pH
    H_deriv = - H * ph_deriv(:,1) * log(10.0)
    block 
        real(r8), dimension(6) :: buffer
        buffer = 0.0_r8
    call check(error, all(is_equal(H_deriv, H_deriv_ref)), .true., " H deriv does not match!")
    if(allocated(error)) return

    buffer = pco2_deriv(:,1)
    call check(error, all(is_equal(buffer, pco2_deriv_ref)), .true., " pCO2 deriv does not match!")
    if(allocated(error)) return

    buffer = fco2_deriv(:,1)
    call check(error, all(is_equal(buffer, fco2_deriv_ref)), .true., " fCO2 deriv does not match!")
    if(allocated(error)) return

    buffer = co2_deriv(:,1)
    call check(error, all(is_equal(buffer, co2_deriv_ref)), .true., " CO2 deriv does not match!")
    if(allocated(error)) return

    buffer = hco3_deriv(:,1)
    call check(error, all(is_equal(buffer, hco3_deriv_ref)), .true., " HCO3 deriv does not match!")
    if(allocated(error)) return

    buffer = co3_deriv(:,1)
    call check(error, all(is_equal(buffer, co3_deriv_ref)), .true., " CO3 deriv does not match!")
    if(allocated(error)) return

    buffer = omegaa_deriv(:,1)
    call check(error, all(is_equal(buffer, omegaa_deriv_ref)), .true., " OmegaA deriv does not match!")
    if(allocated(error)) return

    buffer = omegac_deriv(:,1)
    call check(error, all(is_equal(buffer, omegac_deriv_ref)), .true., " OmegaC deriv does not match!")
    if(allocated(error)) return
    end block

end subroutine test_derivauto

subroutine test_derivnum(error)
   USE mocsy_derivnum
    type(error_type), allocatable, intent(out) :: error 



!  Output variables:
   REAL(kind=rx), DIMENSION(1) :: h, ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis
!  derivative of "vars" Output variables:
   REAL(kind=rx), DIMENSION(1) :: dh_dx, dpco2_dx, dfco2_dx, dco2_dx, dhco3_dx, dco3_dx, dOmegaA_dx, dOmegaC_dx
!  Input variables
   REAL(kind=rx), DIMENSION(1) :: temp, sal, alk, dic, sil, phos, Patm, depth, lat
!  Input options
   CHARACTER(10) :: optCON, optT, optP, optB, optKf, optK1K2
   real(r8), dimension(13) :: h_deriv, pco2_deriv, fco2_deriv, co2_deriv, &
                     hco3_deriv, co3_deriv, omegaa_deriv, omegac_deriv
! Parameter arrays for numerical derivative correctness checking
real(r8), parameter :: H_deriv_num_ref(13) = [ &
    -0.000026065459_r8, 0.000028597328_r8, 0.000029356205_r8, &
    0.000001069695_r8, 0.000000000257_r8, 0.000000000211_r8, &
    0.000000000000_r8, 0.000245081155_r8, 5.785412058112_r8, &
    0.936652884856_r8, 3590.517492074891_r8, 0.000000000000_r8, &
    0.000000000000_r8 ]

real(r8), parameter :: pco2_deriv_num_ref(13) = [ &
    -1214834.992394306464_r8, 1487032.933989916695_r8, 1368207.080659441883_r8, &
    49855.352093665351_r8, 12.826696812900_r8, 8.356155772449_r8, &
    -8994.216689178375_r8, -238075442.011425137520_r8, 231734460834.142639160156_r8, &
    43654657985.432388305664_r8, 167343543999267.250000000000_r8, 0.000000000000_r8, &
    0.000000000000_r8 ]

real(r8), parameter :: fco2_deriv_num_ref(13) = [ &
    -1210607.866089125397_r8, 1481858.670752896927_r8, 1363446.282539268723_r8, &
    49681.875965743471_r8, 12.795486699643_r8, 8.327079785872_r8, &
    -8962.920512919429_r8, -237247037.353988200426_r8, 230928120184.697753906250_r8, &
    43502757723.621490478516_r8, 166761257267041.968750000000_r8, 0.000000000000_r8, &
    0.000000000000_r8 ]

real(r8), parameter :: co2_deriv_num_ref(13) = [ &
    -0.041508848328_r8, 0.050809389671_r8, 0.046749312086_r8, &
    0.001703472703_r8, 0.000000136219_r8, 0.000000227711_r8, &
    0.000000000000_r8, -8.134633489094_r8, 7917.972932260828_r8, &
    1491.605517157524_r8, 5717844.670165921561_r8, 0.000000000000_r8, &
    0.000000000000_r8 ]

real(r8), parameter :: hco3_deriv_num_ref(13) = [ &
    -0.622426414384_r8, 1.575213638168_r8, 0.701007324779_r8, &
    0.025543623826_r8, -0.000000574614_r8, 0.000000934866_r8, &
    0.000000000000_r8, 13.499697422756_r8, -81214.721772032455_r8, &
    22366.672238378316_r8, 85739250.895827606320_r8, 0.000000000000_r8, &
    0.000000000000_r8 ]

real(r8), parameter :: co3_deriv_num_ref(13) = [ &
    0.663934270259_r8, -0.626021939566_r8, -0.747755518413_r8, &
    -0.027247055789_r8, 0.000000438429_r8, -0.000001162563_r8, &
    0.000000000000_r8, -5.365052986271_r8, 73297.007292110007_r8, &
    -23858.242086030546_r8, -91456958.824644848704_r8, 0.000000000000_r8, &
    0.000000000000_r8 ]

real(r8), parameter :: omegaa_deriv_num_ref(13) = [ &
    10253.261127269236_r8, -9667.773912654953_r8, -11547.728341553508_r8, &
    -420.781379221542_r8, 0.016808892531_r8, -0.035806216361_r8, &
    0.000000000000_r8, -82853.516821635349_r8, 1131939394.113049030304_r8, &
    -368447295.314996123314_r8, -1412386922519.985351562500_r8, -4748276.077114485204_r8, &
    0.000000000000_r8 ]

real(r8), parameter :: omegac_deriv_num_ref(13) = [ &
    15855.789291185465_r8, -14950.383509457733_r8, -17857.571859457345_r8, &
    -650.702328145798_r8, 0.014105012177_r8, -0.062941426852_r8, &
    0.000000000000_r8, -128125.860489056955_r8, 1750447228.503089666367_r8, &
    -569772154.134909391403_r8, -2184135287598.605957031250_r8, 0.000000000000_r8, &
    -11355013.399078227580_r8 ]

!  Local variables:
   CHARACTER*3, DIMENSION(13) ::  devar = (/'alk','dic','pho','sil','tem','sal','k0 ','k1 ','k2 ','kb ','kw ','ka ','kc '/)
   INTEGER ::  i
   
  !     derivar = 3-character identifier of input variable with respect to which derivative is requested
  !               possibilities are 'alk', 'dic', 'pho', 'sil', 'tem', or 'sal'
  !

!> Typical options for observations
   optCON  = 'mol/kg'  ! input concentrations are in MOL/KG
   optT    = 'Tinsitu' ! input temperature, variable 'temp' is actually IN SITU temp [°C]
   optP    = 'db'       ! input variable 'depth' is in meters
   optB    = 'l10'
   optK1K2 = 'l'
   optKf   = 'dg'
!> Simple input data (with CONCENTRATION units typical for DATA)
!> (based on observed average surface concentrations from S. Ocean (south of 60°S)--GLODAP and WOA2009)
    temp(1)   = 18.0d0            !Can be "Potential temperature" or "In situ temperature" (see optT below)
    sal(1)    = 35.0d0           !Salinity (practical scale)
    alk(1)    = 2300.0d-6      ! Convert obs. S. Ocean ave surf ALK (umol/kg) to mocsy data units (mol/kg)
    dic(1)    = 2000.0d-6      ! Convert obs. S. Ocean ave surf DIC (umol/kg) to mocsy data units (mol/kg)
    sil(1)    = 0.0d0   ! 60.d-06
    phos(1)   = 0.0d0   !  2.d-06
    sil(1)    = 60.0d-6   ! 60.d-06
    phos(1)   =  2.0d-6   !  2.d-06
    depth(1)  = 0.d0
    Patm(1)   = 1.0d0            !Atmospheric pressure (atm)
    lat(1)    = 0.d0

!  Select input var 'x', choosing set of dy_i/dx to be computed, where y_i are the diff output vars
   call vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,         &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1,                             &  ! INPUT
             optCON, optT, optP, optB=optB, optK1K2=optK1K2, optKf=optKf,                     &
             optGAS='Ppot'    )
   h(1) = 10**(- ph(1))

!    write (*,*) "Variables:"
!    write (*,*) "          h,           ph,         pco2,         fco2,           co2",&
!               "           hco3,           co3,         OmegaA,         OmegaC"
!    write (*,"(9ES15.6)")  h, ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC

!    write (*,*) "Absolute derivatives" 
!    write (*,*) "             dh_dx         dpco2_dx       dfco2_dx         dco2_dx      dhco3_dx       dco3_dx", &
!         "       dOmegaA_dx     dOmegaC_dx"

   do i = 1,13
      call derivnum (dh_dx, dpco2_dx, dfco2_dx, dco2_dx, dhco3_dx,                      &
                      dco3_dx, dOmegaA_dx, dOmegaC_dx,                                   &
                      temp, sal, alk, dic, sil, phos, Patm, depth, lat, 1, devar(i),     &
                      optCON, optT, optP, optB=optB, optK1K2=optK1K2, optKf=optKf          )
      H_deriv(i) = dh_dx(1)
      pco2_deriv(i) = dpco2_dx(1)
      fco2_deriv(i) = dfco2_dx(1)
      co2_deriv(i) = dco2_dx(1)
      hco3_deriv(i) = dhco3_dx(1)
      co3_deriv(i) = dco3_dx(1)
      omegaa_deriv(i) = dOmegaA_dx(1)
      omegac_deriv(i) = dOmegaC_dx(1)
    !   write (*,"(A3,A5,8ES15.6)")  devar(i), "  :  ", dh_dx(1), dpco2_dx(1), dfco2_dx(1), dco2_dx(1), dhco3_dx(1), dco3_dx(1), &
    !       dOmegaA_dx(1), dOmegaC_dx(1)
   end do


   call check(error, all(is_equal(H_deriv, H_deriv_num_ref)), .true., " H deriv does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(pco2_deriv, pco2_deriv_num_ref)), .true., " pCO2 deriv does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(fco2_deriv, fco2_deriv_num_ref)), .true., " fCO2 deriv does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(co2_deriv, co2_deriv_num_ref)), .true., " CO2 deriv does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(hco3_deriv, hco3_deriv_num_ref)), .true., " HCO3 deriv does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(co3_deriv, co3_deriv_num_ref)), .true., " CO3 deriv does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(omegaa_deriv, omegaa_deriv_num_ref)), .true., " OmegaA deriv does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(omegac_deriv, omegac_deriv_num_ref)), .true., " OmegaC deriv does not match!")
   if(allocated(error)) return

end subroutine test_derivnum


subroutine test_errors(error)
   USE mocsy_errors
    type(error_type), allocatable, intent(out) :: error 
    ! Parameter arrays for error correctness checking
! Parameter arrays for error correctness checking
real(r8), parameter :: eh_ref(6) = [ &
    0.000000000381_r8, 0.000000000414_r8, 0.000000000450_r8, &
    0.000000000489_r8, 0.000000000532_r8, 0.000000000579_r8 ]

real(r8), parameter :: epco2_ref(6) = [ &
    19.490615772166_r8, 33.657161894325_r8, 58.343886479187_r8, &
    101.531322754074_r8, 177.382746598312_r8, 311.134506139166_r8 ]

real(r8), parameter :: efco2_ref(6) = [ &
    19.407183076994_r8, 21.869404664938_r8, 24.691263117849_r8, &
    27.931248465059_r8, 31.658153720274_r8, 35.953012514780_r8 ]

real(r8), parameter :: eco2_ref(6) = [ &
    0.000001113102_r8, 0.000001087727_r8, 0.000001064279_r8, &
    0.000001042671_r8, 0.000001022820_r8, 0.000001004647_r8 ]

real(r8), parameter :: ehco3_ref(6) = [ &
    0.000004853454_r8, 0.000004822192_r8, 0.000004791523_r8, &
    0.000004761381_r8, 0.000004731726_r8, 0.000004702533_r8 ]

real(r8), parameter :: eco3_ref(6) = [ &
    0.000003356647_r8, 0.000003347378_r8, 0.000003337425_r8, &
    0.000003326775_r8, 0.000003315435_r8, 0.000003303429_r8 ]

real(r8), parameter :: eOmegaA_ref(6) = [ &
    0.089250247964_r8, 0.072342810004_r8, 0.058883175112_r8, &
    0.048131226039_r8, 0.039512436165_r8, 0.032579510425_r8 ]

real(r8), parameter :: eOmegaC_ref(6) = [ &
    0.141797917915_r8, 0.113526355064_r8, 0.091266038763_r8, &
    0.073677858642_r8, 0.059732550008_r8, 0.048636657333_r8 ]


!  Output variables:
   REAL(kind=rx), DIMENSION(6) :: eh, epco2, efco2, eco2, ehco3, eco3, eOmegaA, eOmegaC
!  Input variables
   REAL(kind=rx), DIMENSION(6) :: temp, sal, alk, dic, sil, phos, Patm, depth, lat
   REAL(kind=rx), DIMENSION(6) :: temp_e, sal_e, ALK_e, DIC_e, sil_e, phos_e
   REAL(kind=rx), DIMENSION(7) :: epK
   
!  Input options
   CHARACTER(10) :: optCON, optT, optP, optB, optKf, optK1K2

!  Local variables:
   INTEGER ::  i

!> Typical options for observations
   optCON  = 'mol/kg'  ! input concentrations are in MOL/KG
   optT    = 'Tinsitu' ! input temperature, variable 'temp' is actually IN SITU temp [°C]
   optP    = 'm'       ! input variable 'depth' is in meters
   optB    = 'u74'
   optK1K2 = 'l'
   optKf   = 'dg'

    ! Input errors
    ALK_e(:) = 2.d-6    ! (2 umol/kg for ALK)
    DIC_e(:) = 2.d-6    ! (2 umol/kg for DIC)
    sal_e(:)  = 0.0d0     ! (psu)
    temp_e(:) = 0.0d0    ! (C)
    phos_e(:) = 0.1d-6
    sil_e(:) = 4.0d-6    
    epK(:) = 0.0

   ! ------------------
   ! 1s test : 1 record 
   ! ------------------
    
    temp(1)   = 18.0d0           ! Can be "Potential temperature" or "In situ temperature" (see optT below)
    sal(1)    = 35.0d0           ! Salinity (practical scale)
    alk(1)    = 2300.d-6         ! Convert obs. S. Ocean ave surf ALK (umol/kg) to mocsy data units (mol/kg)
    dic(1)    = 2000.d-6         ! Convert obs. S. Ocean ave surf DIC (umol/kg) to mocsy data units (mol/kg)
!
    sil(1)    = 60.d-6           ! 60
    phos(1)   = 2.d-6            !  2
!
!   sil(1)    = 0.d0
!   phos(1)   = 0.d0
!
    depth(1)  = 0.d0
    Patm(1)   = 1.0d0            ! Atmospheric pressure (atm)
    lat(1)    = 0.d0

   ! Simple input data (with CONCENTRATION units typical for DATA)
   ! (based on observed average surface concentrations from S. Ocean (south of 60°S)--GLODAP and WOA2009)
   DO i = 1,6
     temp(i)   = 2.0            !Can be "Potential temperature" or "In situ temperature" (see optT below)
     sal(i)    = 35.0           !Salinity (practical scale)
     alk(i)    = 2295.*1.e-6      ! Convert obs. S. Ocean ave surf ALK (umol/kg) to mocsy data units (mol/kg)
     dic(i)    = 2154.*1.e-6      ! Convert obs. S. Ocean ave surf DIC (umol/kg) to mocsy data units (mol/kg)
     sil(i)    = 0.
     phos(i)   = 0.
     depth(i) = real(i-1) * 1000. ! Vary depth from 0 to 5000 db by 1000 db
     Patm(i)   = 1.0            !Atmospheric pressure (atm)
     lat(i)    = 0.
   END DO
   
   ! compute output errors
   call errors(eh, epco2, efco2, eco2, ehco3, eco3, eOmegaA, eOmegaC,          &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, 6,              &  ! INPUT
             temp_e, sal_e, ALK_e, DIC_e, sil_e, phos_e,                       &
             optCON, optT, optP, optB=optB, optK1K2=optK1K2, optKf=optKf   )

    call check(error, all(is_equal(eh, eh_ref)), .true., "eh does not match!")
    if(allocated(error)) return

    call check(error, all(is_equal(epco2, epco2_ref)), .true., "epco2 does not match!")
    if(allocated(error)) return

    call check(error, all(is_equal(efco2, efco2_ref)), .true., "efco2 does not match!")
    if(allocated(error)) return

    call check(error, all(is_equal(eco2, eco2_ref)), .true., "eco2 does not match!")
    if(allocated(error)) return

    call check(error, all(is_equal(ehco3, ehco3_ref)), .true., "ehco3 does not match!")
    if(allocated(error)) return

    call check(error, all(is_equal(eco3, eco3_ref)), .true., "eco3 does not match!")
    if(allocated(error)) return

    call check(error, all(is_equal(eOmegaA, eOmegaA_ref)), .true., "eOmegaA does not match!")
    if(allocated(error)) return

    call check(error, all(is_equal(eOmegaC, eOmegaC_ref)), .true., "eOmegaC does not match!")
    if(allocated(error)) return

end subroutine test_errors

subroutine test_kprime(error)
    use mocsy_gasx
    type(error_type), allocatable, intent(out) :: error
   INTEGER, PARAMETER :: n = 1
  
    real(r8), parameter :: kp_cfc11_ref =    2.1530437923645514E-2_r8
    real(r8), parameter :: kp_cfc12_ref =    5.4134572855389085E-3_r8
    real(r8), parameter :: kp_sf6_ref   =    3.5842379756833637E-4_r8
!  Output variables:
   REAL(kind=r8), DIMENSION(1) :: kp_cfc11, kp_cfc12, kp_sf6
   
!  Input variables
   REAL(kind=rx), DIMENSION(1) :: temp, sal
   
!  Input at standard T and S
!  temp(1)   = 25.0    
!  sal(1)    = 35.0    

   temp(1)   = 4.0    
   sal(1)    = 34.0    

!  Call
   call kprime('cfc11', temp, sal, n, kp_cfc11)
   call kprime('cfc12', temp, sal, n, kp_cfc12)
   call kprime('sf6',   temp, sal, n, kp_sf6)


   call check(error, is_equal(kp_cfc11(1), kp_cfc11_ref), .true., "kp_cfc11 does not match!")
   if(allocated(error)) return

   call check(error, is_equal(kp_cfc12(1), kp_cfc12_ref), .true., "kp_cfc12 does not match!")
   if(allocated(error)) return

   call check(error, is_equal(kp_sf6(1), kp_sf6_ref), .true., "kp_sf6 does not match!")
   if(allocated(error)) return

end subroutine test_kprime

subroutine test_phizero(error)
    use mocsy_gasx
    type(error_type), allocatable, intent(out) :: error
       INTEGER, PARAMETER :: n = 1
   
!  Output variables:
   REAL(kind=r8), DIMENSION(1) :: phi0_cfc11, phi0_cfc12, phi0_sf6, phi0_co2, phi0_n2o
    real(r8), parameter :: ref_phi0_cfc11 =    2.1360647621062642E-002_r8
    real(r8), parameter :: ref_phi0_cfc12 =    5.3702907609728869E-003_r8
    real(r8), parameter :: ref_phi0_sf6   =    3.5607586310647287E-004_r8
    real(r8), parameter :: ref_phi0_co2   =    5.5171775458356262E-002_r8
    real(r8), parameter :: ref_phi0_n2o   =    4.0808900345868551E-002_r8
   
!  Input variables
   REAL(kind=rx), DIMENSION(1) :: temp, sal
   
!  Input at standard T and S
!  temp(1)   = 25.0    
!  sal(1)    = 35.0    

   temp(1)   = 4.0    
   sal(1)    = 34.0    

!  Call
   call phizero('cfc11', temp, sal, n, phi0_cfc11)
   call phizero('cfc12', temp, sal, n, phi0_cfc12)
   call phizero('sf6',   temp, sal, n, phi0_sf6)
   call phizero('co2',   temp, sal, n, phi0_co2)
   call phizero('n2o',   temp, sal, n, phi0_n2o)

   call check(error, is_equal(phi0_cfc11(1), ref_phi0_cfc11), .true., "phi0_cfc11 does not match!")
   if(allocated(error)) return

   call check(error, is_equal(phi0_cfc12(1), ref_phi0_cfc12), .true., "phi0_cfc12 does not match!")
   if(allocated(error)) return

   call check(error, is_equal(phi0_sf6(1), ref_phi0_sf6), .true., "phi0_sf6 does not match!")
   if(allocated(error)) return

   call check(error, is_equal(phi0_co2(1), ref_phi0_co2), .true., "phi0_co2 does not match!")
   if(allocated(error)) return

   call check(error, is_equal(phi0_n2o(1), ref_phi0_n2o), .true., "phi0_n2o does not match!")
   if(allocated(error)) return

end subroutine test_phizero


subroutine test_buffesm(error)
    use mocsy_buffesm
    type(error_type), allocatable, intent(out) :: error

    ! Parameter arrays for output variable correctness checking
real(r8), parameter :: ph_ref(6) = [ &
    8.065774361276_r8, 8.025851077783_r8, 7.985881831551_r8, &
    7.945878996737_r8, 7.905855294833_r8, 7.865823710589_r8 ]

real(r8), parameter :: pco2_ref(6) = [ &
    369.027507401863_r8, 637.882784234534_r8, 1104.109068274922_r8, &
    1913.637800425156_r8, 3321.019885284770_r8, 5770.775448321036_r8 ]

real(r8), parameter :: fco2_ref(6) = [ &
    367.448413463302_r8, 415.986048218755_r8, 471.576625088658_r8, &
    535.308681950995_r8, 608.445696813835_r8, 692.455073657361_r8 ]

real(r8), parameter :: co2_ref(6) = [ &
    0.000021394132_r8, 0.000021030991_r8, 0.000020702172_r8, &
    0.000020405659_r8, 0.000020139600_r8, 0.000019902301_r8 ]

real(r8), parameter :: hco3_ref(6) = [ &
    0.002027862177_r8, 0.002030361926_r8, 0.002032766036_r8, &
    0.002035075500_r8, 0.002037291235_r8, 0.002039414101_r8 ]

real(r8), parameter :: co3_ref(6) = [ &
    0.000104743889_r8, 0.000102607276_r8, 0.000100531980_r8, &
    0.000098519025_r8, 0.000096569346_r8, 0.000094683775_r8 ]

real(r8), parameter :: OmegaA_ref(6) = [ &
    1.577496332906_r8, 1.273001233594_r8, 1.032424280379_r8, &
    0.841530169436_r8, 0.689408747139_r8, 0.567666525488_r8 ]

real(r8), parameter :: OmegaC_ref(6) = [ &
    2.506275274586_r8, 1.997900285257_r8, 1.600618260301_r8, &
    1.288794935890_r8, 1.042978829258_r8, 0.848352925850_r8 ]

real(r8), parameter :: BetaD_ref(6) = [ &
    14.421651930031_r8, 14.485581243375_r8, 14.546107000371_r8, &
    14.603364486957_r8, 14.657466609243_r8, 14.708505263381_r8 ]

real(r8), parameter :: rhoSW_ref(6) = [ &
    1027.971751476259_r8, 1032.629068336757_r8, 1037.187278732428_r8, &
    1041.649029072363_r8, 1046.016875562734_r8, 1050.293287986511_r8 ]

real(r8), parameter :: p_ref(6) = [ &
    0.000000000000_r8, 1000.000000000000_r8, 2000.000000000000_r8, &
    3000.000000000000_r8, 4000.000000000000_r8, 5000.000000000000_r8 ]

real(r8), parameter :: tempis_ref(6) = [ &
    2.000000000000_r8, 2.000000000000_r8, 2.000000000000_r8, &
    2.000000000000_r8, 2.000000000000_r8, 2.000000000000_r8 ]

real(r8), parameter :: gammaDIC_ref(6) = [ &
    0.000149308731_r8, 0.000148673241_r8, 0.000148077249_r8, &
    0.000147518588_r8, 0.000146995414_r8, 0.000146506174_r8 ]

real(r8), parameter :: gammaAlk_ref(6) = [ &
    -0.000166637063_r8, -0.000165743753_r8, -0.000164901357_r8, &
    -0.000164107287_r8, -0.000163359344_r8, -0.000162655671_r8 ]

real(r8), parameter :: betaDIC_ref(6) = [ &
    0.000166637063_r8, 0.000165743753_r8, 0.000164901357_r8, &
    0.000164107287_r8, 0.000163359344_r8, 0.000162655671_r8 ]

real(r8), parameter :: betaAlk_ref(6) = [ &
    -0.000173085147_r8, -0.000172020807_r8, -0.000171012804_r8, &
    -0.000170058534_r8, -0.000169155781_r8, -0.000168302673_r8 ]

real(r8), parameter :: omegaDIC_ref(6) = [ &
    -0.000188515847_r8, -0.000187242983_r8, -0.000186038705_r8, &
    -0.000184899716_r8, -0.000183823198_r8, -0.000182806760_r8 ]

real(r8), parameter :: omegaAlk_ref(6) = [ &
    0.000180052529_r8, 0.000178792203_r8, 0.000177594849_r8, &
    0.000176457819_r8, 0.000175378850_r8, 0.000174356023_r8 ]


!  For vars routine (called below)
!  "vars" Output variables:g
   REAL(kind=rx), DIMENSION(6) :: ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis
!  "vars" Input variables
   INTEGER :: N
   REAL(kind=rx), DIMENSION(6) :: temp, sal, alk, dic, sil, phos, Patm, depth, lat
   REAL(kind=rx), DIMENSION(6)::  gammaDIC, gammaAlk, betaDIC, betaAlk, omegaDIC, omegaAlk, Rf
!  "vars" Input options
   CHARACTER(10) :: optCON, optT, optP, optB, optKf, optK1K2, optGAS
!  CHARACTER(7) :: optGAS

!  Local variables:
   INTEGER :: i


!> Typical options for observations
   optCON  = 'mol/kg'  ! input concentrations are in MOL/KG
   optT    = 'Tinsitu' ! input temperature, variable 'temp' is actually IN SITU temp [°C]
   optP    = 'db'      ! input variable 'depth' is in 'DECIBARS'
   optB    = 'l10'
   optK1K2 = 'l'
   optKf   = 'dg'
   optGAS  = 'Pinsitu'
!> Simple input data (with CONCENTRATION units typical for DATA)
!> (based on observed average surface concentrations from S. Ocean (south of 60°S)--GLODAP and WOA2009)
   DO i = 1,6
     temp(i)   = 2.0            !Can be "Potential temperature" or "In situ temperature" (see optT below)
     sal(i)    = 35.0           !Salinity (practical scale)
     alk(i)    = 2295.*1.e-6      ! Convert obs. S. Ocean ave surf ALK (umol/kg) to mocsy data units (mol/kg)
     dic(i)    = 2154.*1.e-6      ! Convert obs. S. Ocean ave surf DIC (umol/kg) to mocsy data units (mol/kg)
     sil(i)    = 0.
     phos(i)   = 0.
     depth(i) = real(i-1) * 1000. ! Vary depth from 0 to 5000 db by 1000 db
     Patm(i)   = 1.0            !Atmospheric pressure (atm)
     N = i
   END DO

   call vars(ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis,         &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,                             &  ! INPUT
             optCON='mol/kg', optT='Tinsitu', optP='db', optB='l10', optK1K2=optK1K2,         &  ! OPTIONS
             optKf='dg', optGAS=optGAS)                                                          

   call buffesm(gammaDIC, betaDIC, omegaDIC, gammaALK, betaALK, omegaALK, Rf,                    &  ! OUTPUT
             temp, sal, alk, dic, sil, phos, Patm, depth, lat, N,                                 &  ! INPUT
             optCON='mol/kg', optT='Tinsitu', optP='db', optB='l10', optK1K2=optK1K2,             &  ! OPTIONS
             optKf='dg', optGAS=optGAS)     

   call check(error, all(is_equal(ph, ph_ref)), .true., "ph does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(pco2, pco2_ref)), .true., "pCO2 does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(fco2, fco2_ref)), .true., "fCO2 does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(co2, co2_ref)), .true., "CO2 does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(hco3, hco3_ref)), .true., "HCO3 does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(co3, co3_ref)), .true., "CO3 does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(OmegaA, OmegaA_ref)), .true., "OmegaA does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(OmegaC, OmegaC_ref)), .true., "OmegaC does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(BetaD, BetaD_ref)), .true., "BetaD does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(rhoSW, rhoSW_ref)), .true., "rhoSW does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(p, p_ref)), .true., "p does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(tempis, tempis_ref)), .true., "tempis does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(gammaDIC, gammaDIC_ref)), .true., "gammaDIC does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(gammaALK, gammaALK_ref)), .true., "gammaALK does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(betaDIC, betaDIC_ref)), .true., "betaDIC does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(omegaDIC, omegaDIC_ref)), .true., "omegaDIC does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(gammaALK, gammaALK_ref)), .true., "gammaALK does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(betaALK, betaALK_ref)), .true., "betaALK does not match!")
   if(allocated(error)) return

   call check(error, all(is_equal(omegaALK, omegaALK_ref)), .true., "omegaALK does not match!")
   if(allocated(error)) return

end subroutine test_buffesm

end module test_mocsy
