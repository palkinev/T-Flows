!==============================================================================!
  subroutine Residual_Vector(ni, r, b, a, x)
!------------------------------------------------------------------------------!
!   Calculates residual vector {r} = {b} - [A]{x}                              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Matrix_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer           :: ni
  real              :: r(1:)  ! this might be only for inner cells
  real              :: b(1:)  ! this might be only for inner cells
  type(Matrix_Type) :: a
  real              :: x(1:)  ! this may incude buffer cells
!-----------------------------------[Locals]-----------------------------------!
  integer  :: i, j, k, sub
!==============================================================================!

  !----------------!
  !   r = b - Ax   !
  !----------------!
  ! Why not callig this: call exchange(x) ???
  do i = 1, ni
    do j = a % row(i), a % row(i+1) - 1
      k = a % col(j)
      r(i) = b(i) - a % val(j) * x(k)
    end do
  end do

  end subroutine
