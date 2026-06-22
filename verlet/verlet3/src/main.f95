program main
! See subroutines at end of file

! Var definitions, etc
    use kinds, ONLY: wp => dp
    implicit none

    LOGICAL :: OK
    CHARACTER(len=32) :: arg
    CHARACTER(len=32) :: file_name = "atoms.dat"
    integer :: arg_len
    integer :: status
    real (KIND=wp), DIMENSION(:,:), ALLOCATABLE :: x, v, f, fnext, mass
    ! ser, er & der are for parameters to jpca15 function
    real (KIND=wp), DIMENSION(3) :: ser, er, der
    real (KIND=wp), DIMENSION(3) :: d_AB_ser, d_AB_er, d_AB_der
    real (KIND=wp), DIMENSION(3) :: d_AC_ser, d_AC_er, d_AC_der
    real (KIND=wp), DIMENSION(3) :: d_BC_ser, d_BC_er, d_BC_der
    ! for reading in atomic data from file
    integer :: nk
    real :: tau
    real :: sigma
    real :: epsilon
    real :: tmp
    ! interatomic distances
    real (KIND=wp) :: d_AB, d_AC, d_BC
    ! Delta for interatomic distances
    real (KIND=wp), DIMENSION(3)  :: dd_AB, dd_AC, dd_BC
    real (KIND=wp) :: delta = 0.01
    ! For computing jPCA with delta
    real (KIND=wp), DIMENSION(3) :: dd_ser, dd_er, dd_der
    !integer :: num_rows
    integer :: num_atoms
    ! _a & _b are for particles a & b
    real (KIND=wp) :: ax, ay, az, vx, vy, vz
    ! for storing intermediate values
    integer :: i, j, k
    !integer, parameter:: wp = SELECTED_REAL_KIND (p = 13, r = 300)
    integer :: steps = 0
    !real (KIND = wp), DIMENSION(7) :: p_a, p_b ! our two particles


! Vars over
! Program begins
! Process command-line args
    i = 0
    DO
      i = i + 1
      IF (i.gt.command_argument_count()) exit

      CALL get_command_argument(i, arg)
      ! Delta - defaults to 0.01 (see above)
      IF (arg == "-d") THEN
          i = i + 1
          CALL get_command_argument(i, arg)
          print *, "arg -d:", arg
          read( arg, '(f33.32)' ) delta 
          IF (delta == 0) THEN
              print *, "Delta too small!"
              STOP
          END IF
      ! Num steps - defaults to 2000 (see above)
      ELSE IF (arg == "-s") THEN
          i = i + 1
          CALL get_command_argument(i, arg)
          ! print *, "arg -s:", arg
          read( arg, '(I5)' ) steps 
      ! file_name -  defaults to atoms.dat (see above)
      ELSE IF (arg == "-f") THEN
          i = i + 1
          CALL get_command_argument(i, arg)
          file_name = trim(arg)
          INQUIRE (FILE=file_name, EXIST=OK)
          if (.NOT. OK) THEN
              print *, "ERROR!!  File does not exist: ", file_name
              STOP
          END IF
      ELSE
          print *, "Usage: verlet3 [ -h ] [ -f data_file ]  [ -d delta ] [ -s steps ]"
          print *, "DEFAULTS:"
          print '(A, A)', "    data_file - ", file_name
          print '(A, F0.9)', "        delta - ", delta
          print '(A, I0)', "        steps - ", steps
          STOP
      END IF
    END DO

    ! ser = (/0.1, 0.2, 0.3/)
    ! call jpca15(ser, er, der)
    ! print *, er(1), er(2), er(3)
    ! print *, der(1), der(2), der(3)
    ! return

! Read atoms, etc. from data file
    open (UNIT=11, FILE=file_name, STATUS="old", ACTION="read")
    read(unit = 11, FMT=*) nk, tau
    read(unit = 11, FMT=*) sigma, epsilon
    read(unit = 11, FMT=*) num_atoms
    !print *, nk, tau, sigma, epsilon, num_atoms
    print *, "Number of atoms:", num_atoms

    ! Let the command-line `steps` override nk, if it's set
    if (steps > 0) then
        nk = steps
    end if

    print *, "Will use delta: ", delta
    print *, "Will use file: ", file_name
    print *, "Will use num steps: ", nk

! Allocate arrays for position, velocity, force & mass
    ! Position
    allocate(x(num_atoms,3))
    ! Velocity
    allocate(v(num_atoms,3))
    ! Force
    allocate(f(num_atoms,3))
    ! Force
    allocate(fnext(num_atoms,3))
    ! mass
    allocate(mass(num_atoms,1))

    ! read in info for particles
    do i = 1, num_atoms, 1
        read(unit = 11, FMT=*) mass(i, 1), ax, ay, az, vx, vy, vz
        x(i,:) = (/ ax, ay, az /)
        v(i,:) = (/ vx, vy, vz /)
        print *, 'Particle', i, ': '
        print *, '  Mass:', mass(i, 1)
        print *, '  Starting Position:', x(i, :)
        print *, '  Initial Velocity: ', v(i, :)
    end do
    close(unit = 11) 

! Initial Force for particle a
    do a = 1, num_atoms, 1
        do d = 1, 3, 1 ! dimensions
            f(a, d) = lj(epsilon, sigma, x(1,:), x(2,:), 1)
        end do
    end do
    f(1,1) = lj(epsilon, sigma, x(1,:), x(2,:), 1)
    f(1,2) = lj(epsilon, sigma, x(1,:), x(2,:), 2)
    f(1,3) = lj(epsilon, sigma, x(1,:), x(2,:), 3)
    f(2,1) = lj(epsilon, sigma, x(2,:), x(1,:), 1)
    f(2,2) = lj(epsilon, sigma, x(2,:), x(1,:), 2)
    f(2,3) = lj(epsilon, sigma, x(2,:), x(1,:), 3)
    print *, 'initial f(1, :) ', f(1, :)
    print *, 'initial f(2, :) ', f(2, :)

  ! Compute 
  do i = 1, nk, 1
    ! euclidean distances
    !d_AB = SQRT( (x(1, 1) - x(2, 1))**2 +  (x(1, 2) - x(2, 2))**2 +  (x(1, 3) - x(2, 3))**2 ) 
    call eudist( x(1, :) , x(2, :), d_AB)
    call eudist( x(1, :) , x(3, :), d_BC)
    call eudist( x(2, :) , x(3, :), d_AC)
    ! Initial values of potentials
    ser = (/d_AB, d_AC, d_BC/)
    call jpca15(ser, er, der)
    ! Numerical derivative of each inter-atomic distance, in each direction
    ! E.g., AB x, y, z
    ! Call subroutine eudist_with_delta(p1, p2, delta, dim, dist)
    call eudist_with_delta(x(1, :), x(2, :), delta, 1, dd_AB(1))
    call eudist_with_delta(x(1, :), x(2, :), delta, 2, dd_AB(2))
    call eudist_with_delta(x(1, :), x(2, :), delta, 3, dd_AB(3))
    call eudist_with_delta(x(1, :), x(3, :), delta, 1, dd_AC(1))
    call eudist_with_delta(x(1, :), x(3, :), delta, 2, dd_AC(2))
    call eudist_with_delta(x(1, :), x(3, :), delta, 3, dd_AC(3))
    call eudist_with_delta(x(2, :), x(3, :), delta, 1, dd_BC(1))
    call eudist_with_delta(x(2, :), x(3, :), delta, 2, dd_BC(2))
    call eudist_with_delta(x(2, :), x(3, :), delta, 3, dd_BC(3))
    ! print *, "eudist_with_delta - dd_BC(3)", dd_BC(3)

!    print *, "er: ", er
!    print *, "der: ", der
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
! replace with jpca15
!    fnext = lj(epsilon, sigma, x)

    v(1,1) = v(1,1) + tau/(2*mass(1,1)) * (f(1,1) + fnext(1,1) )
    v(1,2) = v(1,2) + tau/(2*mass(1,1)) * (f(1,2) + fnext(1,2) )
    v(1,3) = v(1,3) + tau/(2*mass(1,1)) * (f(1,3) + fnext(1,3) )
    v(2,1) = v(2,1) + tau/(2*mass(2,1)) * (f(2,1) + fnext(2,1) )
    v(2,2) = v(2,2) + tau/(2*mass(2,1)) * (f(2,2) + fnext(2,2) )
    v(2,3) = v(2,3) + tau/(2*mass(2,1)) * (f(2,3) + fnext(2,3) )

  end do

  print *, 'p1: ', x(1, :)
  print *, 'p2: ', x(2, :)
  print *, 'p3: ', x(3, :)

! Cleanup
    deallocate(x)
    deallocate(v)
    deallocate(f)
    deallocate(fnext)

end program main

subroutine eudist(p1, p2, dist)
    use kinds, ONLY: wp => dp
    ! Two points
    real (KIND=wp), DIMENSION(3) :: p1, p2
    ! Distance to return
    real (KIND=wp) :: dist
    dist = SQRT( (p1(1) - p2(1))**2 +  (p1(2) - p2(2))**2 +  (p1(3) - p2(3))**2 ) 
end subroutine eudist

subroutine eudist_with_delta(p1, p2, delta, dim, dist)
    use kinds, ONLY: wp => dp
    ! Two points
    real (KIND=wp), DIMENSION(3) :: p1, p2
    ! delta to add to dimension, and distance to return
    real (KIND=wp) :: delta, dist
    ! Dimension which to add delta
    integer :: dim
    p1(dim) = p1(dim) + delta
    p2(dim) = p2(dim) + delta
    call eudist(p1, p2, dist)
    ! Not divide by delta to get "derivative"
    dist = dist / delta
end subroutine eudist_with_delta
