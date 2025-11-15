module lennard_jones
  USE kinds, ONLY: wp => dp
  implicit none
  
  PUBLIC :: lj

  CONTAINS

  function dist(x, y, z) result(r)
    real :: x, y, z, r
    r = sqrt(x**2 + y**2 + z**2)
  end function dist

  function v_prime_lj(epsilon, sigma, r) result(v_prime)
    real, intent(in) :: epsilon, sigma, r
    real :: v_prime
    v_prime = 4 * epsilon * (-12 * sigma**12 / r**13 + 6 * sigma**6 / r**7 )
  end function v_prime_lj

  function lj(epsilon, sigma, particle_a, particle_b) result(f_a)
    real, intent(in) :: epsilon, sigma
    real (KIND=wp), INTENT(IN) :: particle_a, particle_b
    real :: f_a
    !f_a = - (dist(particle_a, particle_b))
  end function lj
 
end module lennard_jones
