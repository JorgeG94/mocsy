!> \file singledouble.F90
!! \BRIEF 
!> Module that defines single and double precision - used by all other modules
MODULE mocsy_singledouble

IMPLICIT NONE ; PRIVATE

PUBLIC rx, r8, wp
public :: sgle


#if USE_PRECISION == 2
  INTEGER, PARAMETER :: rx = KIND(1.0d0)
#else
  INTEGER, PARAMETER :: rx = KIND(1.0)
#endif

  INTEGER, PARAMETER :: r8 = KIND(1.0d0)
  integer, parameter :: r4 = kind(1.0)
  INTEGER, PARAMETER :: wp = KIND(1.0d0)

  contains 

    elemental FUNCTION SGLE(x) RESULT(sgle_result)
    !> Convert to appropriate precision based on USE_PRECISION
    !> In double precision mode (USE_PRECISION==2): returns input unchanged
    !> In single precision mode: converts to single precision
    !> \param x value to convert
    !> \return precision-adjusted value
    REAL(rx), INTENT(IN) :: x
    REAL(rx) :: sgle_result

#if USE_PRECISION == 2
    ! Double precision mode - return as-is
    sgle_result = x
#else
    ! Single precision mode - convert to single precision
    sgle_result = REAL(x, KIND=rx)
#endif

  END FUNCTION SGLE

END MODULE mocsy_singledouble
