program main
  use kinds, only: wp => dp
  use mymodule
  implicit none
  integer :: i
  real (kind=wp), DIMENSION(:), ALLOCATABLE :: x, y
  i = 3
  ALLOCATE ( x(i), y(i) )
  x = (/ 1.0, 2.0, 5.0 /)
  CALL myexpsubroutine(x, i, y)
  print *, y
  DEALLOCATE (x, y)
end program main
