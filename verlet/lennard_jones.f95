module lennard_jones
  USE kinds, ONLY: wp => dp
  implicit none
  
  PUBLIC :: lj, v_prime_lj

  CONTAINS

  real function v_prime_lj(epsilon, sigma, r_ab) result(v_prime)
    real, intent(in) :: epsilon, sigma
    real, intent(in) :: r_ab
    !real, intent(out) :: v_prime
    v_prime = 4.0_wp * epsilon * (-12.0_wp * sigma**12.0_wp / r_ab**13.0_wp + 6.0_wp* sigma**6.0_wp / r_ab**7.0_wp)
  end function v_prime_lj

  logical function realeq(a, b, tol)
    ! function to check "equaility" of reals
    implicit none
    real, intent(in) :: a, b, tol
    !print *, 'realeq: comparing ', abs(a - b), tol
    realeq = abs(a - b) < tol
  end function

  logical function is0(a)
    ! function to check if a number is 0 - or pretty close
    implicit none
    real, intent(in) :: a
    real :: tol
    tol = 1E-35
    is0 = realeq(a, 1E-36, tol)
  end function

  real function r_ab(v) result(r)
    real (KIND=wp), dimension(3), INTENT(IN) :: v
    r = sqrt( v(1)**2 + v(2)**2 + v(3)**2 )
  end function r_ab

  function distance(a, b) result(d)
    implicit none
    real (KIND=wp), dimension(3), INTENT(IN) :: a, b
    real (KIND=wp), dimension(3) :: d
    d = (/ a(1) - b(1), a(2) - b(2), a(3) - b(3)   /)
  end function

  function lj(epsilon, sigma, particles) result(f)
    real, intent(in) :: epsilon, sigma
    real (KIND=wp), dimension(2,3), INTENT(IN) :: particles
    real (KIND=wp), dimension(2,3) :: f
    !real :: f_a
    real (KIND=wp), dimension(3) :: dist_ab, dist_ba
    real :: r
    real (KIND=wp) :: vplj_a, vplj_b
    dist_ab = distance( particles(1,:), particles(2,:) )
    dist_ba = distance( particles(2,:), particles(1,:) )
    vplj_a = v_prime_lj(epsilon, sigma, r_ab(dist_ab) )
    vplj_b = v_prime_lj(epsilon, sigma, r_ab(dist_ba) )
    r = r_ab(dist_ab) ! either direction will do
    f(1,:) = - (/ ((dist_ab(1) / r) * vplj_a), ((dist_ab(2) / r) * vplj_a), ((dist_ab(3) / r) * vplj_a)/)
    f(2,:) = - (/ ((dist_ba(1) / r) * vplj_a), ((dist_ba(2) / r) * vplj_a), ((dist_ba(3) / r) * vplj_a)/)
   !print *,  'lj: ', f_a, axes(axis), r_ab, epsilon, sigma
  end function lj
 
end module lennard_jones
