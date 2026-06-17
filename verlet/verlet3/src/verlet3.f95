program main

    use kinds, ONLY: wp => dp
    implicit none

    CHARACTER(len=32) :: arg
    CHARACTER(len=32) :: file_name = "atoms.dat"
    integer :: arg_len
    integer :: status
    real (KIND=wp), DIMENSION(:,:), ALLOCATABLE :: x, v, f, fnext, mass
    real (KIND=wp), DIMENSION(3) :: ser, er, der
    real (KIND=wp), DIMENSION(3) :: dx_AB_ser, dx_AB_er, dx_AB_der
    real (KIND=wp), DIMENSION(3) :: dy_AB_ser, dy_AB_er, dy_AB_der
    real (KIND=wp), DIMENSION(3) :: dz_AB_ser, dz_AB_er, dz_AB_der
    real (KIND=wp), DIMENSION(3) :: dx_AC_ser, dx_AC_er, dx_AC_der
    real (KIND=wp), DIMENSION(3) :: dy_AC_ser, dy_AC_er, dy_AC_der
    real (KIND=wp), DIMENSION(3) :: dz_AC_ser, dz_AC_er, dz_AC_der
    real (KIND=wp), DIMENSION(3) :: dx_BC_ser, dx_BC_er, dx_BC_der
    real (KIND=wp), DIMENSION(3) :: dy_BC_ser, dy_BC_er, dy_BC_der
    real (KIND=wp), DIMENSION(3) :: dz_BC_ser, dz_BC_er, dz_BC_der
    integer :: nk
    real :: tau
    real :: sigma
    real :: epsilon
    real :: tmp
    ! interatomic distances
    real (KIND=wp) :: d_AB, d_AC, d_BC
    ! Delta for interatomic distances
    real (KIND=wp) :: dx_AB, dy_AB, dz_AB
    real (KIND=wp) :: dx_AC, dy_AC, dz_AC
    real (KIND=wp) :: dx_BC, dy_BC, dz_BC
    real (KIND=wp) :: delta = 0.000001
    ! For computing jPCA with delta
    real (KIND=wp), DIMENSION(3) :: dd_ser, dd_er, dd_der
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
      i = i + 1
      IF (i.gt.command_argument_count()) exit

      CALL get_command_argument(i, arg)
      IF (arg == "-d") THEN
          i = i + 1
          CALL get_command_argument(i, arg)
          print *, "arg -d:", arg
          read( arg, '(f33.32)' ) delta 
          IF (delta == 0) THEN
              print *, "Delta too small!"
              STOP
          END IF
      ELSE IF (arg == "-f") THEN
          i = i + 1
          CALL get_command_argument(i, arg)
          file_name = trim(arg)
      ELSE
          print *, "Usage: verlet3 [ [-f data_file] | atoms.dat ]  [ [-d delta |", delta, "]"
          STOP
      END IF
    END DO
    print *, "Will use delta: ", delta
    print *, "Will use file: ", file_name

    ! ser = (/0.1, 0.2, 0.3/)
    ! call jpca15(ser, er, der)
    ! print *, er(1), er(2), er(3)
    ! print *, der(1), der(2), der(3)
    ! return

    open (UNIT=11, FILE=file_name, STATUS="old", ACTION="read")
    read(unit = 11, FMT=*) nk, tau
    read(unit = 11, FMT=*) sigma, epsilon
    read(unit = 11, FMT=*) num_atoms
    !print *, nk, tau, sigma, epsilon, num_atoms
    print *, "Number of atoms:", num_atoms

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
        print *, 'Particle', i, ': '
        print *, '  Starting Position:', x(i, :)
        print *, '  Initial Velocity: ', v(i, :)
    end do
    close(unit = 11) 
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
    print *, 'A', x(1, :)
    print *, 'B', x(2, :)
    print *, 'C', x(3, :)

  do i = 1, steps, 1
    ! euclidean distances
    d_AB = SQRT( (x(1, 1) - x(2, 1))**2 +  (x(1, 2) - x(2, 2))**2 +  (x(1, 3) - x(2, 3))**2 ) 
    d_AC = SQRT( (x(1, 1) - x(3, 1))**2 +  (x(1, 2) - x(3, 2))**2 +  (x(1, 3) - x(3, 3))**2 ) 
    d_BC = SQRT( (x(2, 1) - x(3, 1))**2 +  (x(2, 2) - x(3, 2))**2 +  (x(2, 3) - x(3, 3))**2 ) 
    ! initial values of potentials
    ser = (/d_AB, d_AC, d_BC/)
    call jpca15(ser, er, der)
    ! numerical derivative of each inter-atomic distance, in each direction
    ! AB x, y, z
    dx_AB = SQRT( (x(1, 1) - x(2, 1) + delta)**2 + (x(1, 2) - x(2, 2))**2 +  (x(1, 3) - x(2, 3))**2 ) 
    dy_AB = SQRT( (x(1, 1) - x(2, 1))**2 + (x(1, 2) - x(2, 2) + delta)**2 +  (x(1, 3) - x(2, 3))**2 ) 
    dz_AB = SQRT( (x(1, 1) - x(2, 1))**2 + (x(1, 2) - x(2, 2))**2 +  (x(1, 3) - x(2, 3) + delta)**2 ) 
    ! AC x, y, z
    dx_AC = SQRT( (x(1, 1) - x(3, 1) + delta)**2 + (x(1, 2) - x(3, 2))**2 +  (x(1, 3) - x(3, 3))**2 ) 
    dy_AC = SQRT( (x(1, 1) - x(3, 1))**2 + (x(1, 2) - x(3, 2) + delta)**2 +  (x(1, 3) - x(3, 3))**2 ) 
    dz_AC = SQRT( (x(1, 1) - x(3, 1))**2 + (x(1, 2) - x(3, 2))**2 +  (x(1, 3) - x(3, 3) + delta)**2 ) 
    ! BC x, y, z
    dx_BC = SQRT( (x(2, 1) - x(3, 1) + delta)**2 + (x(2, 2) - x(3, 2))**2 +  (x(2, 3) - x(3, 3))**2 ) 
    dy_BC = SQRT( (x(2, 1) - x(3, 1))**2 + (x(2, 2) - x(3, 2) + delta)**2 +  (x(2, 3) - x(3, 3))**2 ) 
    dz_BC = SQRT( (x(2, 1) - x(3, 1))**2 + (x(2, 2) - x(3, 2))**2 +  (x(2, 3) - x(3, 3) + delta)**2 ) 

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

    ! if ( MOD(i,100) .EQ. 0 ) then
    !   print *, i
    !   print *, 'p1 x: ', x(1, :)
    !   print *, 'p1 v: ', v(1, :)
    !   print *, 'p1 f: ', f(1, :)
! !      print *, 'p2: ', x(2, :)
    ! end if

  end do

  print *, 'p1: ', x(1, :)
  print *, 'p2: ', x(2, :)
  print *, 'p3: ', x(3, :)

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
