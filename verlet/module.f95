module mymodule
  USE kinds, ONLY: wp => dp
  implicit none
  
  PUBLIC :: myexpsubroutine

  CONTAINS

  subroutine myexpsubroutine(x, i, y)
    integer, INTENT(in) :: i
    real (KIND=wp), DIMENSION(:), INTENT(IN) :: x
    real (KIND=wp), DIMENSION(:), INTENT(OUT) :: Y
    y = x**i
  end subroutine myexpsubroutine
 
end module mymodule
