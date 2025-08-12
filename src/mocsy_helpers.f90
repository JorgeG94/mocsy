!! General helpers
module mocsy_test_helpers
  !! simple reusable helpers for random things
   use mocsy_singledouble, only: rx
   implicit none

   private

   real(rx), parameter :: tol_dp = 1.0e-12_rx
   public :: is_equal

   interface is_equal
      module procedure is_equal_dp
   end interface is_equal

contains

   elemental function is_equal_dp(a, b, tol) result(res)
      real(rx), intent(in) :: a, b
      real(rx), intent(in), optional :: tol
      real(rx) :: working_tol
      logical :: res

      if(present(tol)) then 
        working_tol = tol
      else 
        working_tol = tol_dp
      end if

      res = abs(a - b) < working_tol
   end function is_equal_dp

end module mocsy_test_helpers
