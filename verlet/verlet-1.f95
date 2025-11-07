! Write a Fortran program that implements the Verlet algorithm with $k$ ranging
! from 1 to 10 and $\tau$ = 0.2 s for one particle of mass 1 kg in 3D space subject
! to a constant force expressed by components
! $$f^{(a, x)}$$ = 0 kg m s-2
! $$f^{(a, y)}$$ = 0.1 kg m s-2
! $$f^{(a, z)}$$ = 0 kg m s-2.

program verlet
  IMPLICIT NONE

  integer, parameter:: wp = SELECTED_REAL_KIND (p = 13, r = 300)

  real (KIND = wp):: tau = 0.2
  real (KIND = wp):: mass = 1.0
  ! Force functions?
  ! real (kind = wp):: fx = 0, fy = 0.1, fz = 0
  real (kind = wp):: f = 1
  integer:: i
  integer, parameter:: steps = 10
  ! lists of stuff?
  real (KIND = wp), dimension(steps+1):: a_x, a_y, a_z
  real (KIND = wp), dimension(steps+1):: v_x, v_y, v_z
  real (KIND = wp), dimension(steps+1):: f_x, f_y, f_z

  a_x(1) = 0.0_wp
  a_y(1) = 0.0_wp
  a_z(1) = 0.0_wp

  v_x(1) = 0.0_wp
  v_y(1) = 0.0_wp
  v_z(1) = 0.0_wp

  f_x(1) = 0.0_wp
  f_y(1) = 0.1_wp
  f_z(1) = 0.0_wp

  do i = 2, steps, 1
    ! \alpha_{k+1}
    a_x(i) = a_x(i-1) + tau*v_x(i-1) + tau**2 * f_x(i-1) / (2*mass)
    a_y(i) = a_y(i-1) + tau*v_y(i-1) + tau**2 * f_y(i-1) / (2*mass)
    a_z(i) = a_z(i-1) + tau*v_z(i-1) + tau**2 * f_z(i-1) / (2*mass)

    ! evaluate f, but f is constant for now
    f_x(i) = f_x(i-1)
    f_y(i) = f_y(i-1)
    f_z(i) = f_z(i-1)

    ! v_{k+1}
    v_x(i) = v_x(i-1) + tau/2*mass * (f_x(i-1) + f_x(i) )
    v_y(i) = v_y(i-1) + tau/2*mass * (f_y(i-1) + f_y(i) )
    v_z(i) = v_z(i-1) + tau/2*mass * (f_z(i-1) + f_z(i) )
    print *, a_x(i), a_y(i), a_z(i)
  end do
end program verlet
