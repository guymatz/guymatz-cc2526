module lennard_jones
  USE kinds, ONLY: wp => dp
  implicit none
  
  PUBLIC :: lj

  CONTAINS

!  function v_prime_lj(epsilon, sigma, r_ab) result(v_prime)
!    real, intent(in) :: epsilon, sigma
!    real, intent(in) :: r_ab
!    real, intent(out) :: v_prime
!    !v_prime = 4 * epsilon * (-12 * sigma**12 / r_ab**13 + 6 * sigma**6 / r_ab**7 )
!  end function v_prime_lj

  function lj(epsilon, sigma, particle_a, particle_b, axis) result(f_a)
    integer :: axis
    real, intent(in) :: epsilon, sigma
    real (KIND=wp), dimension(3), INTENT(IN) :: particle_a, particle_b
    real :: f_a
    real :: axes(3)
    real :: x_ab, y_ab, z_ab
    real :: r_ab
    x_ab = particle_a(1) - particle_b(1)
    y_ab = particle_a(2) - particle_b(2)
    z_ab = particle_a(3) - particle_b(3)
    axes = (/ x_ab, y_ab, z_ab /)
    r_ab = sqrt( x_ab**2 + y_ab**2 + z_ab**2 )
    !f_a = - ( axes(axis) / r_ab ) * v_prime_lj(epsilon, sigma, r_ab)

  end function lj
 
end module lennard_jones
