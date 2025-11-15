program main
    
    use lennard_jones
    use kinds, ONLY: wp => dp
    implicit none

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
    real (KIND = wp), dimension(steps+1):: a_x_a, a_y_a, a_z_a, a_x_b, a_y_b, a_z_b
    real (KIND = wp), dimension(steps+1):: v_x_a, v_y_a, v_z_a, v_x_b, v_y_b, v_z_b
    real (KIND = wp), dimension(steps+1):: f_x_a, f_y_a, f_z_a, f_x_b, f_y_b, f_z_b

    open (UNIT=11, FILE="verlet-2.dat", STATUS="old", ACTION="read")
    ! num_rows = count_lines(11)
    read(unit = 11, FMT=*) nk, tau
    read(unit = 11, FMT=*) sigma, epsilon
    read(unit = 11, FMT=*) num_atoms
    !print *, nk, tau, sigma, epsilon, num_atoms
    print *, num_atoms

    read(unit = 11, FMT=*) m, x, y, z, vx, vy, vz
    p_a = (/ m, x, y, z, vx, vy, vz /)
    read(unit = 11, FMT=*) m, x, y, z, vx, vy, vz
    p_b = (/ m, x, y, z, vx, vy, vz /)
    close(unit = 11) 
!    print *, 'Particle A:'
!    print *, p_a
!    print *, 'Particle B:'
!    print *, p_b
!    print *


    ! Initial position for particle a
    a_x_a(1) = p_a(2)
    a_y_a(1) = p_a(3)
    a_z_a(1) = p_a(4)
    ! Initial Velocity for particle a
    v_x_a(1) = p_a(5)
    v_y_a(1) = p_a(6)
    v_z_a(1) = p_a(7)

    ! Initial Force for particle a
    f_x_a(1) = 0.0_wp
    f_y_a(1) = 0.0_wp
    f_z_a(1) = 0.0_wp

    ! same for particle b
    ! Initial position for particle a
    a_x_b(1) = p_b(2)
    a_y_b(1) = p_b(3)
    a_z_b(1) = p_b(4)
    ! Initial Velocity for particle a
    v_x_b(1) = p_b(5)
    v_y_b(1) = p_b(6)
    v_z_b(1) = p_b(7)

    ! Initial Force for particle a
    f_x_b(1) = 0.0_wp
    f_y_b(1) = 0.0_wp
    f_z_b(1) = 0.0_wp


    print *, '1'
    print *, 'a', a_x_a(1), a_y_a(1), a_z_a(1)
    print *, 'b', a_x_b(1), a_y_b(1), a_z_b(1)

  do i = 2, steps, 1
    print *, i
    print *, 'a', a_x_a(i), a_y_a(i), a_z_a(i)
    print *, 'b', a_x_b(i), a_y_b(i), a_z_b(i)
    ! \alpha_{k+1}
    a_x_a(i) = a_x_a(i-1) + tau*v_x_a(i-1) + tau**2 * f_x_a(i-1) / (2*p_a(1))
    a_y_a(i) = a_y_a(i-1) + tau*v_x_a(i-1) + tau**2 * f_x_a(i-1) / (2*p_a(1))
    a_z_a(i) = a_z_a(i-1) + tau*v_x_a(i-1) + tau**2 * f_x_a(i-1) / (2*p_a(1))
    ! particle b
    a_x_b(i) = a_x_b(i-1) + tau*v_x_b(i-1) + tau**2 * f_x_b(i-1) / (2*p_b(1))
    a_y_b(i) = a_y_b(i-1) + tau*v_x_b(i-1) + tau**2 * f_x_b(i-1) / (2*p_b(1))
    a_z_b(i) = a_z_b(i-1) + tau*v_x_b(i-1) + tau**2 * f_x_b(i-1) / (2*p_b(1))

    ! evaluate f on particle a
    f_x_a(i) = lj(epsilon, sigma, p_a, p_b, 1)
    f_x_a(i) = lj(epsilon, sigma, p_a, p_b, 2)
    f_x_a(i) = lj(epsilon, sigma, p_a, p_b, 3)

    ! evaluate f on particle b
    f_x_b(i) = lj(epsilon, sigma, p_b, p_a, 1)
    f_x_b(i) = lj(epsilon, sigma, p_b, p_a, 2)
    f_x_b(i) = lj(epsilon, sigma, p_b, p_a, 3)

    ! v_{k+1}
    v_x_a(i) = v_x_a(i-1) + tau/2*p_a(1) * (f_x_a(i-1) + f_x_a(i) )
    v_x_a(i) = v_x_a(i-1) + tau/2*p_a(1) * (f_x_a(i-1) + f_x_a(i) )
    v_x_a(i) = v_x_a(i-1) + tau/2*p_a(1) * (f_x_a(i-1) + f_x_a(i) )
  end do
!  print *, a_x_a(i), a_y_a(i), a_z_a(i)

! Write a Fortran program that implements the Verlet algorithm with $k$ ranging
! from 1 to 10 and $\tau$ = 0.2 s for one particle of mass 1 kg in 3D space subject
! to a constant force expressed by components
! $$f^{(a, x)}$$ = 0 kg m s-2
! $$f^{(a, y)}$$ = 0.1 kg m s-2
! $$f^{(a, z)}$$ = 0 kg m s-2.



end program main
