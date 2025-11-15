program main
    
    use lennard_jones
    use kinds, ONLY: wp => dp

    integer :: nk
    real :: tau
    real :: sigma
    real :: epsilon
    integer :: num_rows
    integer :: num_atoms
    real :: m, x, y, z, vx, vy, vz
    integer :: i, j, k
    !integer, parameter:: wp = SELECTED_REAL_KIND (p = 13, r = 300)
    integer, parameter:: steps = 10
    real (KIND = wp), DIMENSION(7) :: p_a, p_b ! our two particles
  ! lists of stuff?
    real (KIND = wp), dimension(steps+1):: a_x, a_y, a_z
    real (KIND = wp), dimension(steps+1):: v_x, v_y, v_z
    real (KIND = wp), dimension(steps+1):: f_x, f_y, f_z

    open (UNIT=11, FILE="verlet-2.dat", STATUS="old", ACTION="read")
    ! num_rows = count_lines(11)
    read(unit = 11, FMT=*) nk, tau
    read(unit = 11, FMT=*) sigma, epsilon
    read(unit = 11, FMT=*) num_atoms
    print *, nk, tau, sigma, epsilon, num_atoms

    read(unit = 11, FMT=*) m, x, y, z, vx, vy, vz
    p_a = (/ m, x, y, z, vx, vy, vz /)
    read(unit = 11, FMT=*) m, x, y, z, vx, vy, vz
    p_b = (/ m, x, y, z, vx, vy, vz /)
    close(unit = 11) 
    print *, 'Particle A:'
    print *, p_a
    print *, 'Particle B:'
    print *, p_b


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
    !print *, a_x(i), a_y(i), a_z(i)
    ! \alpha_{k+1}
    a_x(i) = a_x(i-1) + tau*v_x(i-1) + tau**2 * f_x(i-1) / (2*mass)
    a_y(i) = a_y(i-1) + tau*v_y(i-1) + tau**2 * f_y(i-1) / (2*mass)
    a_z(i) = a_z(i-1) + tau*v_z(i-1) + tau**2 * f_z(i-1) / (2*mass)

    ! evaluate f, but f is constant for now
    f_x(i) = lj(epsilon, sigma, p_a(1), p_b(1))
    f_y(i) = lj(epsilon, sigma, p_a(2), p_b(2))
    f_z(i) = lj(epsilon, sigma, p_a(3), p_b(3))

    ! v_{k+1}
    v_x(i) = v_x(i-1) + tau/2*mass * (f_x(i-1) + f_x(i) )
    v_y(i) = v_y(i-1) + tau/2*mass * (f_y(i-1) + f_y(i) )
    v_z(i) = v_z(i-1) + tau/2*mass * (f_z(i-1) + f_z(i) )
  end do
  print *, a_x(i), a_y(i), a_z(i)

! Write a Fortran program that implements the Verlet algorithm with $k$ ranging
! from 1 to 10 and $\tau$ = 0.2 s for one particle of mass 1 kg in 3D space subject
! to a constant force expressed by components
! $$f^{(a, x)}$$ = 0 kg m s-2
! $$f^{(a, y)}$$ = 0.1 kg m s-2
! $$f^{(a, z)}$$ = 0 kg m s-2.



end program main
