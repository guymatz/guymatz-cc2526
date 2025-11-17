program main
    
    use lennard_jones
    use kinds, ONLY: wp => dp
    implicit none

    real (KIND=wp), DIMENSION(:,:), ALLOCATABLE :: x, v, f, fnext
    integer :: nk
    real :: tau
    real :: sigma
    real :: epsilon
    real :: tmp
    !integer :: num_rows
    integer :: num_atoms
    real (KIND=wp) :: mass, ax, ay, az, vx, vy, vz
    ! for storing intermediate values
    integer :: i!, j, k
    !integer, parameter:: wp = SELECTED_REAL_KIND (p = 13, r = 300)
    integer, parameter:: steps = 10
    !real (KIND = wp), DIMENSION(7) :: p_a, p_b ! our two particles


    allocate(x(2,3))
    allocate(v(2,3))
    allocate(f(2,3))
    allocate(fnext(2,3))
    open (UNIT=11, FILE="verlet-2.dat", STATUS="old", ACTION="read")
    ! num_rows = count_lines(11)
    read(unit = 11, FMT=*) nk, tau
    read(unit = 11, FMT=*) sigma, epsilon
    read(unit = 11, FMT=*) num_atoms
    !print *, nk, tau, sigma, epsilon, num_atoms
    print *, num_atoms

    ! read in info for particles a & b
    read(unit = 11, FMT=*) mass, ax, ay, az, vx, vy, vz
    x(1,:) = (/ 0.0, 1.1, 2.2 /)
    !x(1,:) = (/ ax, ay, az /)
    v(1,:) = (/ vx, vy, vz /)
    read(unit = 11, FMT=*) mass, ax, ay, az, vx, vy, vz
    x(2,:) = (/ ax, ay, az /)
    v(2,:) = (/ vx, vy, vz /)
    close(unit = 11) 
!    print *, 'Particle A:'
!    print *, p_a
!    print *, 'Particle B:'
!    print *, p_b
!    print *
    ! Initial Force for particle a
    f(1,:) = (/ 0.0_wp, 0.0_wp, 0.0_wp /)
    f(2,:) = (/ 0.0_wp, 0.0_wp, 0.0_wp /)


    print *, '1'
    print *, 'a', x(1, :)
    print *, 'b', x(2, :)

  do i = 2, steps, 1
    ! \alpha_{k+1}
    x(1,1) = x(1,1) + tau*v(1,1) + tau**2 * f(1,1) / (2*x(1,1))
    x(1,2) = x(1,2) + tau*v(1,2) + tau**2 * f(1,2) / (2*x(1,2))
    x(1,3) = x(1,3) + tau*v(1,3) + tau**2 * f(1,3) / (2*x(1,3))
    ! particle b
    x(2,1) = x(2,1) + tau*v(2,1) + tau**2 * f(2,1) / (2*x(2,1))
    x(2,2) = x(2,2) + tau*v(2,2) + tau**2 * f(2,2) / (2*x(2,2))
    x(2,3) = x(2,3) + tau*v(2,3) + tau**2 * f(2,3) / (2*x(2,3))
    print *, 'x 1', x(1, 1),  x(1, 2),  x(1, 3)
    print *, 'x 2', x(2, 1),  x(2, 2),  x(2, 3)

    ! evaluate f on particle a
    !print *, 'lj x: ', epsilon, sigma, p_a, p_b, 1
    fnext(1,1) = lj(epsilon, sigma, x(1,:), x(2,:), 1)
    fnext(1,2) = lj(epsilon, sigma, x(1,:), x(2,:), 2)
    fnext(1,3) = lj(epsilon, sigma, x(1,:), x(2,:), 3)
    ! evaluate fnext on particle b
    fnext(2,1) = lj(epsilon, sigma, x(2,:), x(1,:), 1)
    fnext(2,2) = lj(epsilon, sigma, x(2,:), x(1,:), 2)
    fnext(2,3) = lj(epsilon, sigma, x(2,:), x(1,:), 3)

    ! v_{k+1}
    v(1,1) = v(1,1) + tau/2*x(1,1) * (f(1,1) + fnext(1,1) )
    v(1,2) = v(1,2) + tau/2*x(1,2) * (f(1,2) + fnext(1,2) )
    v(1,3) = v(1,3) + tau/2*x(1,3) * (f(1,3) + fnext(1,3) )
    v(2,1) = v(2,1) + tau/2*x(2,1) * (f(2,1) + fnext(2,1) )
    v(2,2) = v(2,2) + tau/2*x(2,2) * (f(2,2) + fnext(2,2) )
    v(2,3) = v(2,3) + tau/2*x(2,3) * (f(2,3) + fnext(2,3) )

    f(1,1) = fnext(1,1)
    f(1,2) = fnext(1,2)
    f(1,3) = fnext(1,3)
    ! evaluate fnext on particle b
    f(2,1) = fnext(2,1)
    f(2,2) = fnext(2,2)
    f(2,3) = fnext(2,3)
  end do
  print *, 'p1: ', x(1, :)
  print *, 'p2: ', x(2, :)
deallocate(x)
deallocate(v)
deallocate(f)
deallocate(fnext)

! Write a Fortran program that implements the Verlet algorithm with $k$ ranging
! from 1 to 10 and $\tau$ = 0.2 s for one particle of mass 1 kg in 3D space subject
! to a constant force expressed by components
! $$f^{(a, x)}$$ = 0 kg m s-2
! $$f^{(a, y)}$$ = 0.1 kg m s-2
! $$f^{(a, z)}$$ = 0 kg m s-2.



end program main
