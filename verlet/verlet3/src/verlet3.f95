program main
    use jpca15
!    use lennard_jones
    use kinds, ONLY: wp => dp
    implicit none

    CHARACTER(len=32) :: arg
    real (KIND=wp), DIMENSION(:,:), ALLOCATABLE :: x, v, f, fnext, mass
    real (KIND=wp), DIMENSION(3) :: ser, er, der
    integer :: nk
    real :: tau
    real :: sigma
    real :: epsilon
    real :: tmp
    !integer :: num_rows
    integer :: num_atoms
    ! _a & _b are for particles a & b
    real (KIND=wp) :: ax, ay, az, vx, vy, vz
    ! for storing intermediate values
    integer :: i, j, k
    !integer, parameter:: wp = SELECTED_REAL_KIND (p = 13, r = 300)
    integer, parameter:: steps = 2000
    !real (KIND = wp), DIMENSION(7) :: p_a, p_b ! our two particles

! Vars over
! Program begins
! Process command-line args

    i = 0
    DO
      CALL get_command_argument(i, arg)
      IF (LEN_TRIM(arg) == 0) EXIT
      WRITE (*, *) TRIM(arg)
      i = i + 1
    END DO

    STOP
!    ser = (/0.1, 0.2, 0.3/)
!    call jpca15(ser, er, der)
!    print *, er(1), er(2), er(3)
!    print *, der(1), der(2), der(3)
!    return

    open (UNIT=11, FILE="atoms.dat", STATUS="old", ACTION="read")
    read(unit = 11, FMT=*) nk, tau
    read(unit = 11, FMT=*) sigma, epsilon
    read(unit = 11, FMT=*) num_atoms
    !print *, nk, tau, sigma, epsilon, num_atoms
    print *, num_atoms

    allocate(x(num_atoms,3))
    allocate(v(num_atoms,3))
    allocate(f(num_atoms,3))
    allocate(fnext(num_atoms,3))
    allocate(mass(num_atoms,1))

    ! read in info for particles a & b
    do i = 1, num_atoms, 1
        read(unit = 11, FMT=*) mass(i, 1), ax, ay, az, vx, vy, vz
        x(i,:) = (/ ax, ay, az /)
        v(i,:) = (/ vx, vy, vz /)
    end do
    close(unit = 11) 
!    print *, 'Particle A:'
!    print *, p_a
!    print *, 'Particle B:'
!    print *, p_b
!    print *
    ! Initial Force for particle a
!    f(1,1) = lj(epsilon, sigma, x(1,:), x(2,:), 1)
!    f(1,2) = lj(epsilon, sigma, x(1,:), x(2,:), 2)
!    f(1,3) = lj(epsilon, sigma, x(1,:), x(2,:), 3)
!    f(2,1) = lj(epsilon, sigma, x(2,:), x(1,:), 1)
!    f(2,2) = lj(epsilon, sigma, x(2,:), x(1,:), 2)
!    f(2,3) = lj(epsilon, sigma, x(2,:), x(1,:), 3)
!    print *, 'initial f(1, :) ', f(1, :)
!    print *, 'initial f(2, :) ', f(2, :)

    print *, '1'
    print *, 'a', x(1, :)
    print *, 'b', x(2, :)

  do i = 1, steps, 1
    f = lj(epsilon, sigma, x)
    ! \alpha_{k+1}
    x(1,1) = x(1,1) + tau*v(1,1) + tau**2 * f(1,1) / (2*mass(1, 1))
    x(1,2) = x(1,2) + tau*v(1,2) + tau**2 * f(1,2) / (2*mass(1, 1))
    x(1,3) = x(1,3) + tau*v(1,3) + tau**2 * f(1,3) / (2*mass(1, 1))
    ! particle b
    x(2,1) = x(2,1) + tau*v(2,1) + tau**2 * f(2,1) / (2*mass(2, 1))
    x(2,2) = x(2,2) + tau*v(2,2) + tau**2 * f(2,2) / (2*mass(2, 1))
    x(2,3) = x(2,3) + tau*v(2,3) + tau**2 * f(2,3) / (2*mass(2, 1))
    !print *, 'x 1', i, x(1, 1),  x(1, 2),  x(1, 3)
    !print *, 'x 2', i, x(2, 1),  x(2, 2),  x(2, 3)

    ! evaluate fnext on particle a
    fnext = lj(epsilon, sigma, x)

    ! v_{k+1}
    v(1,1) = v(1,1) + tau/(2*mass(1,1)) * (f(1,1) + fnext(1,1) )
    v(1,2) = v(1,2) + tau/(2*mass(1,1)) * (f(1,2) + fnext(1,2) )
    v(1,3) = v(1,3) + tau/(2*mass(1,1)) * (f(1,3) + fnext(1,3) )
    v(2,1) = v(2,1) + tau/(2*mass(2,1)) * (f(2,1) + fnext(2,1) )
    v(2,2) = v(2,2) + tau/(2*mass(2,1)) * (f(2,2) + fnext(2,2) )
    v(2,3) = v(2,3) + tau/(2*mass(2,1)) * (f(2,3) + fnext(2,3) )

!    f(1,1) = fnext(1,1)
!    f(1,2) = fnext(1,2)
!    f(1,3) = fnext(1,3)
!    ! evaluate fnext on particle b
!    f(2,1) = fnext(2,1)
!    f(2,2) = fnext(2,2)
!    f(2,3) = fnext(2,3)
    if ( MOD(i,100) .EQ. 0 ) then
      print *, i
      print *, 'p1 x: ', x(1, :)
      print *, 'p1 v: ', v(1, :)
      print *, 'p1 f: ', f(1, :)
!      print *, 'p2: ', x(2, :)
    end if
  end do
  print *, 'p1: ', x(1, :)
!  print *, 'p2: ', x(2, :)
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
