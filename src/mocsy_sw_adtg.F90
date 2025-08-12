!> \file sw_adtg.F90
!! \BRIEF 
!> Module with sw_adtg function - compute adiabatic temp. gradient from S,T,P
MODULE mocsy_sw_adtg

USE mocsy_singledouble, only : rx, r8, wp
USE Dual_Num_Auto_Diff

IMPLICIT NONE ; PRIVATE

PUBLIC sw_adtg, sw_adtg_DNAD

CONTAINS
!>  Function to calculate adiabatic temperature gradient as per UNESCO 1983 routines.
FUNCTION sw_adtg(s, t, p)
  !     ==================================================================
  !     Calculates adiabatic temperature gradient as per UNESCO 1983 routines.
  !     Armin Koehl akoehl@ucsd.edu
  !     ==================================================================
  !> salinity [psu (PSU-78)]
  REAL(kind=r8) :: s
  !> temperature [degree C (IPTS-68)]
  REAL(kind=r8) :: t
  !> pressure [db]
  REAL(kind=r8) :: p
  REAL(kind=r8) :: sw_adtg
  REAL(kind=r8), PARAMETER :: sref = 35.0d0
  REAL(kind=r8), PARAMETER :: a0 = 3.5803d-5
  REAL(kind=r8), PARAMETER :: a1 = 8.5258d-6
  REAL(kind=r8), PARAMETER :: a2 = -6.836d-8
  REAL(kind=r8), PARAMETER :: a3 = 6.6228d-10
  REAL(kind=r8), PARAMETER :: b0 = 1.8932d-6
  REAL(kind=r8), PARAMETER :: b1 = -4.2393d-8
  REAL(kind=r8), PARAMETER :: c0 = 1.8741d-8
  REAL(kind=r8), PARAMETER :: c1 = -6.7795d-10
  REAL(kind=r8), PARAMETER :: c2 = 8.733d-12
  REAL(kind=r8), PARAMETER :: c3 = -5.4481d-14
  REAL(kind=r8), PARAMETER :: d0 = -1.1351d-10
  REAL(kind=r8), PARAMETER :: d1 = 2.7759d-12
  REAL(kind=r8), PARAMETER :: e0 = -4.6206d-13
  REAL(kind=r8), PARAMETER :: e1 = 1.8676d-14
  REAL(kind=r8), PARAMETER :: e2 = -2.1687d-16
  

  
  ! UNESCO 1983 adiabatic temperature gradient coefficients
  sw_adtg = a0 + (a1 + (a2 + a3*t)*t)*t &
          + (b0 + b1*t)*(s - sref) &
          + ((c0 + (c1 + (c2 + c3*t)*t)*t) + (d0 + d1*t)*(s - sref))*p &
          + (e0 + (e1 + e2*t)*t)*p*p
          
END FUNCTION sw_adtg

!>  Function to calculate adiabatic temperature gradient as per UNESCO 1983 routines.
!! and derivative with respect to temperature and salinity
FUNCTION sw_adtg_DNAD  (s,t,p)

  !     ==================================================================
  !     Calculates adiabatic temperature gradient as per UNESCO 1983 routines.
  !     Armin Koehl akoehl@ucsd.edu
  !     ==================================================================

  !> salinity [psu (PSU-78)]
  TYPE(DUAL_NUM) :: s
  !> temperature [degree C (IPTS-68)]
  TYPE(DUAL_NUM) :: t
  !> pressure [db]
  TYPE(DUAL_NUM) :: p

  REAL(kind=r8), PARAMETER :: sref = 35.0d0
  REAL(kind=r8), PARAMETER :: a0 = 3.5803d-5
  REAL(kind=r8), PARAMETER :: a1 = 8.5258d-6
  REAL(kind=r8), PARAMETER :: a2 = -6.836d-8
  REAL(kind=r8), PARAMETER :: a3 = 6.6228d-10
  REAL(kind=r8), PARAMETER :: b0 = 1.8932d-6
  REAL(kind=r8), PARAMETER :: b1 = -4.2393d-8
  REAL(kind=r8), PARAMETER :: c0 = 1.8741d-8
  REAL(kind=r8), PARAMETER :: c1 = -6.7795d-10
  REAL(kind=r8), PARAMETER :: c2 = 8.733d-12
  REAL(kind=r8), PARAMETER :: c3 = -5.4481d-14
  REAL(kind=r8), PARAMETER :: d0 = -1.1351d-10
  REAL(kind=r8), PARAMETER :: d1 = 2.7759d-12
  REAL(kind=r8), PARAMETER :: e0 = -4.6206d-13
  REAL(kind=r8), PARAMETER :: e1 = 1.8676d-14
  REAL(kind=r8), PARAMETER :: e2 = -2.1687d-16
  

  TYPE(DUAL_NUM) :: sw_adtg_DNAD

  sw_adtg_DNAD =  a0 + (a1 + (a2 + a3*T)*T)*T &
       + (b0 + b1*T)*(S-sref) &
       + ( (c0 + (c1 + (c2 + c3*T)*T)*T) + (d0 + d1*T)*(S-sref) )*P &
       + (  e0 + (e1 + e2*T)*T )*P*P

END FUNCTION sw_adtg_DNAD
END MODULE mocsy_sw_adtg
