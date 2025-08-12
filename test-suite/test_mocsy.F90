module test_mocsy
  use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
  use, intrinsic :: iso_fortran_env, only: error_unit
   USE mocsy_singledouble
   USE mocsy_constants
   USE mocsy_vars
   USE mocsy_derivauto
   use mocsy_test_helpers, only: is_equal
   use mocsy_reference_values !! all reference values are here!
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


!  Computed variables:
   REAL(kind=r8), DIMENSION(1) :: k0_co2, k0_n2o

!  Input variables
   REAL(kind=rx), DIMENSION(1) :: temp, sal


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
