program main
! See subroutines at end of file
    ! eudist_with_delta(p1, p2, delta, dim, dist)
    ! eudist(p1, p2, dist)

! Var definitions, etc
    use kinds, ONLY: wp => dp
    implicit none
    ! for looping by atoms & dimension
    integer :: atom_num, dim

    LOGICAL :: OK
    CHARACTER(len=32) :: arg
    CHARACTER(len=32) :: file_name = "atoms.dat"
    integer :: arg_len
    integer :: status
    real(KIND=wp), DIMENSION(:, :), ALLOCATABLE :: x, v, f, fnext, mass
    ! ser, er & der are for parameters to jpca15 function
    ! INPUT
    !   ser: a vector with the three interatomic distances (AB, AC, and BC)
    ! OUTPUT
    !   er: potential energy (in eV)
    !   der: vector of the derivatives of the potential with
    !        respect to the interatomic distances (AB, AC, and BC) (in bohr)
    real(KIND=wp), DIMENSION(3) :: d_AB_ser, d_AB_er, d_AB_der
    real(KIND=wp), DIMENSION(3) :: d_AC_ser, d_AC_er, d_AC_der
    real(KIND=wp), DIMENSION(3) :: d_BC_ser, d_BC_er, d_BC_der
    ! Force vectors (in the x-direction only)
    real(KIND=wp), DIMENSION(3) :: f_AB, f_AC, f_BC, force
    ! for reading in atomic data from file
    integer :: nk
    real :: sigma, epsilon
    real :: tau
    real :: tmp
    ! interatomic distances
    real(KIND=wp) :: d_AB, d_AC, d_BC
    ! Delta for interatomic distances
    real(KIND=wp), DIMENSION(3) :: dd_AB, dd_AC, dd_BC
    real(KIND=wp) :: delta = 0.1, cli_delta = 0
    ! For computing jPCA with delta
    real(KIND=wp), DIMENSION(3) :: dd_ser, dd_er, dd_der
    !integer :: num_rows
    integer :: num_atoms
    ! locations and velocites for atom
    real(KIND=wp) :: ax, ay, az, vx, vy, vz
    ! for storing intermediate values & looping
    integer :: i, j, k
    !integer, parameter:: wp = SELECTED_REAL_KIND (p = 13, r = 300)
    integer :: steps = 0
    !real (KIND = wp), DIMENSION(7) :: p_a, p_b ! our two particles

! Process command-line args
    i = 0
    DO
        i = i + 1
        IF (i .gt. command_argument_count()) exit

        CALL get_command_argument(i, arg)
        ! Delta - defaults to 0.01 (see above)
        IF (arg == "-d") THEN
            i = i + 1
            CALL get_command_argument(i, arg)
            print *, "arg -d:", arg
            read (arg, '(f33.32)') cli_delta
            IF (cli_delta == 0) THEN
                print *, "Delta too small!"
                STOP
            END IF
            ! Num steps - defaults to 2000 (see above)
        ELSE IF (arg == "-s") THEN
            i = i + 1
            CALL get_command_argument(i, arg)
            ! print *, "arg -s:", arg
            read (arg, '(I5)') steps
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

! Read atoms, etc. from data file
    open (UNIT=11, FILE=file_name, STATUS="old", ACTION="read")
    read (unit=11, FMT=*) nk, tau
    read (unit=11, FMT=*) delta
    read (unit=11, FMT=*) num_atoms
    !print *, nk, tau, sigma, epsilon, num_atoms
    print *, "Number of atoms:", num_atoms

    ! Let the command-line `steps` override nk, if it's set
    if (steps > 0) then
        nk = steps
    end if
    ! same for delta
    if (cli_delta > 0) then
        delta = cli_delta
    end if

    print *, "Will use delta: ", delta
    print *, "Will use file: ", file_name
    print *, "Will use num steps: ", nk

    ! Allocate arrays for position, velocity, force & mass
    ! Position
    allocate (x(num_atoms, 3))
    ! Velocity
    allocate (v(num_atoms, 3))
    ! Force
    allocate (f(num_atoms, 3))
    ! Force
    allocate (fnext(num_atoms, 3))
    ! mass
    allocate (mass(num_atoms, 1))

    ! read in info for particles
    do i = 1, num_atoms, 1
        read (unit=11, FMT=*) mass(i, 1), ax, ay, az, vx, vy, vz
        x(i, :) = (/ax, ay, az/)
        v(i, :) = (/vx, vy, vz/)
        print *, 'Particle', i, ': '
        print *, '  Mass:', mass(i, 1)
        print *, '  Starting Position:', x(i, :)
        print *, '  Initial Velocity: ', v(i, :)
    end do
    close (unit=11)

! Initial Force for particles
    ! jpca15 function:
    ! INPUT
    !   ser: a vector with the three interatomic distances (AB, AC, and BC)
    ! OUTPUT
    !   er: potential energy (in eV)
    !   der: vector of the derivatives of the potential with
    !        respect to the interatomic distances (AB, AC, and BC) (in bohr)
    call compute_force(x, delta, fnext)

    !d_AB_ser =
    do atom_num = 1, num_atoms, 1
        do dim = 1, 3, 1 ! dimensions

        end do
    end do
    print *, 'initial f_AB ', f_AB
    print *, 'initial f_AC ', f_AC
    print *, 'initial f_BC ', f_BC

    STOP 1

    ! Iterate!
    do k = 1, nk, 1
        ! Calculate x^{(a)}_{k+1}
        do atom_num = 1, num_atoms, 1
            do dim = 1, 3, 1 ! dimensions
                x(atom_num, dim) = x(atom_num, dim) + tau * v(atom_num, dim) + tau**2 * f(atom_num, dim) / (2 * mass(atom_num, dim))
                v(atom_num, dim) = v(atom_num, dim) + tau / (2 * mass(atom_num, dim)) * (f(atom_num, dim) + fnext(atom_num, dim))
            end do
        end do
        ! Calculate f^{(a, x)}_{k+1}
        ! Calculate v^{(a, x)}_{k+1}

    end do

    print *, 'p1: ', x(1, :)
    print *, 'p2: ', x(2, :)
    print *, 'p3: ', x(3, :)

! Cleanup
    deallocate (x)
    deallocate (v)
    deallocate (f)
    deallocate (fnext)

end program main

!function calc_x(point,) result()
!end function calc_x

subroutine eudist(p1, p2, dist)
    use kinds, ONLY: wp => dp
    implicit none
    ! Two points
    real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2
    ! Distance to return
    real(KIND=wp), intent(out) :: dist
    dist = SQRT((p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2)
end subroutine eudist

subroutine get_ser(p1, p2, p3, ser)
    use kinds, ONLY: wp => dp
    implicit none
    ! Three points
    real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2, p3
    real(KIND=wp), DIMENSION(3), intent(out) :: ser

    call eudist(p1, p2, ser(1))
    call eudist(p1, p3, ser(2))
    call eudist(p2, p2, ser(3))
end subroutine get_ser

subroutine get_delta_ser(p1, p2, p3, delta, ser)
    use kinds, ONLY: wp => dp
    implicit none
    ! points
    real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2, p3
    real(KIND=wp), intent(in) :: delta
    ! return n x d (3) vector
    real(KIND=wp), DIMENSION(3, 3), intent(out) :: ser
    ! for looping
    integer :: dim

    do dim = 1, 3, 1
        call eudist_with_delta(p1, p2, delta, dim, ser(1, dim))
        call eudist_with_delta(p1, p3, delta, dim, ser(2, dim))
        call eudist_with_delta(p2, p2, delta, dim, ser(3, dim))
    end do
end subroutine get_delta_ser

subroutine eudist_with_delta(p1, p2, delta, dim, dist)
    use kinds, ONLY: wp => dp
    implicit none
    ! Two points
    real(KIND=wp), DIMENSION(3) :: p1, p2
    ! delta to add to dimension, and distance to return
    real(KIND=wp), intent(in) :: delta
    ! Dimension which to add delta
    integer, intent(in) :: dim
    ! distance to return
    real(KIND=wp), intent(out) :: dist
    p1(dim) = p1(dim) + delta
    !!!!!   NOT Adding delta to second atom.  Only the first (above)
    !p2(dim) = p2(dim) + delta
    !call eudist(p1, p2, dist)
    dist = SQRT((p1(dim) - p2(dim))**2)

    ! print *, "p1: ", p1
    ! print *, "p2: ", p2
    ! print *, "dim: ", dim
    ! print *, "dist: ", dist

    ! Not divide by delta to get "derivative"
    ! dist = dist / delta
end subroutine eudist_with_delta

!    call compute_force(x, delta f_AC)
! points = x,y,z coords of each point
subroutine compute_force(points, delta, force)
    use kinds, ONLY: wp => dp
    implicit none
    integer :: n
    ! particles
    real(KIND=wp), DIMENSION(3, 3), intent(in) :: points
    real(KIND=wp), intent(in) :: delta
    ! distances between particles
    real(KIND=wp) :: d_AB, d_AC, d_BC
    ! Vector of forces betweem each particle
    real(KIND=wp), DIMENSION(3, 3), intent(out) :: force
    !
    real(KIND=wp), DIMENSION(3) :: ser, der
    real(KIND=wp) :: er
    real(KIND=wp), DIMENSION(3) :: delta_ser, delta_der
    real(KIND=wp) :: delta_er
    ! n*(n-1)/2 gets number of pairs of particles
    real(KIND=wp), DIMENSION(3) :: d_
    ! n*(n-1)/2 gets number of pairs of particles in each direction
    real(KIND=wp), DIMENSION(3, 3) :: delta_d_
    ! delta to add to dimension, and distance to return
    ! for looping
    integer :: dim, atom_i, atom_j

    ! get der first for actual location
    ! First we get interatomic distances
    call eudist(points(1, :), points(2, :), d_(1))
    call eudist(points(1, :), points(3, :), d_(2))
    call eudist(points(2, :), points(3, :), d_(3))
    ! now ser for location + delta

    ser = (/d_(1), d_(2), d_(3)/)
    call jpca15(ser, er, der)

    do dim = 1, 3, 1
        call eudist_with_delta(points(1, :), points(2, :), delta, dim, delta_d_(1, dim))
        call eudist_with_delta(points(1, :), points(3, :), delta, dim, delta_d_(2, dim))
        call eudist_with_delta(points(2, :), points(3, :), delta, dim, delta_d_(3, dim))
    end do

    do dim = 1, 3, 1
        do atom_i = 1, 3, 1
            delta_ser = (/delta_d_(1, dim), delta_d_(2, dim), delta_d_(3, dim)/)
            call jpca15(delta_ser, delta_er, delta_der)
            force(atom_i, dim) = (delta_der(dim) - der(dim)) / delta
        end do
    end do

    print *, "d_(1): ", d_(1)
    print *, "d_(2): ", d_(2)
    print *, "d_(3): ", d_(3)

    print *, "delta_d_(1): ", delta_d_(1, :)
    print *, "delta_d_(2): ", delta_d_(2, :)
    print *, "delta_d_(3): ", delta_d_(3, :)

    print *, "force_(1): ", force(1, :)
    print *, "force_(2): ", force(2, :)
    print *, "force_(3): ", force(3, :)
    stop 3

    ! ser = (/d_AB, d_AC, d_BC/)
    ! call jpca15(ser, er, der)
    ! call eudist_with_delta(p1, p3, delta, 1, d_AC)
    ! call eudist_with_delta(p2, p3, delta, 1, d_BC)
    ! delta_ser = (/d_AB, d_AC, d_BC/)
    ! call jpca15(delta_ser, delta_er, delta_der)

end subroutine compute_force
