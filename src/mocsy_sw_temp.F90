!> \file sw_temp.F90
!! \BRIEF 
!> Module with sw_temp function - compute in-situ T from potential T
MODULE mocsy_sw_temp

USE mocsy_singledouble, only : rx, r8, wp
USE mocsy_sw_ptmp, only : compute_sea_water_potential_temperature
USE Dual_Num_Auto_Diff

IMPLICIT NONE ; PRIVATE

PUBLIC compute_sea_water_insitu_temperature

interface compute_sea_water_insitu_temperature
  module procedure :: compute_sea_water_insitu_temperature
  module procedure :: compute_sea_water_insitu_temperature_DNAD
end interface compute_sea_water_insitu_temperature

CONTAINS
!> Function to compute in-situ temperature [C] from potential temperature [C]
FUNCTION compute_sea_water_insitu_temperature( s, t, p, pr ) result(sea_water_temperature)
  !     =============================================================
  !     SW_TEMP
  !     Computes in-situ temperature [C] from potential temperature [C]
  !     Routine available in seawater.f (used for MIT GCM)
  !     Downloaded seawater.f (on 17 April 2009) from
  !     http://ecco2.jpl.nasa.gov/data1/beaufort/MITgcm/bin/
  !     =============================================================

  !     REFERENCES:
  !     Fofonoff, P. and Millard, R.C. Jr
  !     Unesco 1983. Algorithms for computation of fundamental properties of
  !     seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
  !     Eqn.(31) p.39

  !     Bryden, H. 1973.
  !     "New Polynomials for thermal expansion, adiabatic temperature gradient
  !     and potential temperature of sea water."
  !     DEEP-SEA RES., 1973, Vol20,401-408.
  !     =============================================================

  !     Simple modifications: J. C. Orr, 16 April 2009
  !     - combined fortran code from MITgcm site & simplification in
  !       CSIRO code (matlab equivalent) from Phil Morgan

#if USE_PRECISION == 2
#   define SGLE(x)    (x)
#else
#   define SGLE(x)    REAL(x)
#endif

  !     Input arguments:
  !     -----------------------------------------------
  !     s  = salinity              [psu      (PSS-78) ]
  !     t  = potential temperature [degree C (IPTS-68)]
  !     p  = pressure              [db]
  !     pr = reference pressure    [db]

  !> salinity [psu (PSS-78)]
  REAL(kind=rx) ::   s
  !> potential temperature [degree C (IPTS-68)]
  REAL(kind=rx) ::   t
  !> pressure [db]
  REAL(kind=rx) ::   p
  !> reference pressure [db]
  REAL(kind=rx) ::   pr

  REAL(kind=r8) ::  ds, dt, dp, dpr
  REAL(kind=r8) :: dsw_temp

  REAL(kind=rx) ::   sea_water_temperature
! EXTERNAL sw_ptmp
! REAL(kind=r8) ::   sw_ptmp

  ds = DBLE(s)
  dt = DBLE(t)
  dp = DBLE(p)
  dpr = DBLE(pr)

  !    Simple solution
  !    (see https://svn.mpl.ird.fr/us191/oceano/tags/V0/lib/matlab/seawater/sw_temp.m)
  !    Carry out inverse calculation by swapping P_ref (pr) and Pressure (p)
  !    in routine that is normally used to compute potential temp from temp
  dsw_temp = compute_sea_water_potential_temperature(ds, dt, dpr, dp)
  sea_water_temperature = SGLE(dsw_temp)

  !    The above simplification works extremely well (compared to Table in 1983 report)
  !    whereas the sw_temp routine from MIT GCM site does not seem to work right

  RETURN
END FUNCTION compute_sea_water_insitu_temperature


!> Function to compute in-situ temperature [C] from potential temperature [C]
!! and derivative with respect to potential temperature and salinity
FUNCTION compute_sea_water_insitu_temperature_DNAD( s, t, p, pr ) result(sea_water_temperature)
  !     It is similar to subroutine 'sw_temp' above except that it also computes
  !     partial derivative of insitu temperature
  !     with respect to potential temperature and salinity.
  !     =============================================================
  !     SW_TEMP
  !     Computes in-situ temperature [C] from potential temperature [C]
  !     Routine available in seawater.f (used for MIT GCM)
  !     Downloaded seawater.f (on 17 April 2009) from
  !     http://ecco2.jpl.nasa.gov/data1/beaufort/MITgcm/bin/
  !     =============================================================

  !     REFERENCES:
  !     Fofonoff, P. and Millard, R.C. Jr
  !     Unesco 1983. Algorithms for computation of fundamental properties of
  !     seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
  !     Eqn.(31) p.39

  !     Bryden, H. 1973.
  !     "New Polynomials for thermal expansion, adiabatic temperature gradient
  !     and potential temperature of sea water."
  !     DEEP-SEA RES., 1973, Vol20,401-408.
  !     =============================================================

  !     Simple modifications: J. C. Orr, 16 April 2009
  !     - combined fortran code from MITgcm site & simplification in
  !       CSIRO code (matlab equivalent) from Phil Morgan

  !     Input arguments:
  !     -----------------------------------------------
  !     s  = salinity              [psu      (PSS-78) ]
  !     t  = potential temperature [degree C (IPTS-68)]
  !     p  = pressure              [db]
  !     pr = reference pressure    [db]

  !> salinity [psu (PSS-78)]
  TYPE(DUAL_NUM) ::   s
  !> potential temperature [degree C (IPTS-68)]
  TYPE(DUAL_NUM) ::   t
  !> pressure [db]
  TYPE(DUAL_NUM) ::   p
  !> reference pressure [db]
  TYPE(DUAL_NUM) ::   pr

  TYPE(DUAL_NUM) ::   sea_water_temperature


  !    Simple solution
  !    (see https://svn.mpl.ird.fr/us191/oceano/tags/V0/lib/matlab/seawater/sw_temp.m)
  !    Carry out inverse calculation by swapping P_ref (pr) and Pressure (p)
  !    in routine that is normally used to compute potential temp from temp
  sea_water_temperature = compute_sea_water_potential_temperature(s, t, pr, p)

  !    The above simplification works extremely well (compared to Table in 1983 report)
  !    whereas the sw_temp routine from MIT GCM site does not seem to work right

  RETURN
END FUNCTION compute_sea_water_insitu_temperature_DNAD
END MODULE mocsy_sw_temp
