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
    real(KIND=wp), DIMENSION(3) :: f_ABx, f_ACx, f_BCx, force
    ! for reading in atomic data from file
    integer :: nk
    real :: sigma, epsilon
    real :: tau
    real :: tmp
    ! interatomic distances
    real(KIND=wp) :: d_AB, d_AC, d_BC
    ! Delta for interatomic distances
    real(KIND=wp), DIMENSION(3) :: dd_AB, dd_AC, dd_BC
    real(KIND=wp) :: delta = 0.1, cli_delta=0
    ! For computing jPCA with delta
    real(KIND=wp), DIMENSION(3) :: dd_ser, dd_er, dd_der
    !integer :: num_rows
    integer :: num_atoms
    ! _a & _b are for particles a & b
    real(KIND=wp) :: ax, ay, az, vx, vy, vz
    ! for storing intermediate values
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
    call compute_force( x(1, :), x(2, :), x(3, :), delta, f_ABx)
    call compute_force( x(1, :), x(2, :), x(3, :), delta, f_ACx)
    call compute_force( x(1, :), x(2, :), x(3, :), delta, f_BCx)

    !d_AB_ser = 
    do atom_num = 1, num_atoms, 1
        do dim = 1, 3, 1 ! dimensions

        end do
    end do
    print *, 'initial f_ABx ', f_ABx
    print *, 'initial f_ACx ', f_ACx
    print *, 'initial f_BCx ', f_BCx

    ! Compute
    do k = 1, nk, 1
        ! print *, "eudist_with_delta - dd_BC(3)", dd_BC(3)
        do atom_num = 1, num_atoms, 1
            do dim = 1, 3, 1 ! dimensions
                x(atom_num, dim) = x(atom_num, dim) + tau * v(atom_num, dim) + tau**2 * f(atom_num, dim) / (2 * mass(atom_num, dim))
                v(atom_num, dim) = v(atom_num, dim) + tau / (2 * mass(atom_num, dim)) * (f(atom_num, dim) + fnext(atom_num, dim))
            end do
        end do

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

subroutine eudist(p1, p2, dist)
    use kinds, ONLY: wp => dp
    implicit none
    ! Two points
    real(KIND=wp), DIMENSION(3) :: p1, p2
    ! Distance to return
    real(KIND=wp) :: dist
    dist = SQRT((p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2)
end subroutine eudist

subroutine get_ser(p1, p2, p3, ser)
    use kinds, ONLY: wp => dp
    implicit none
    ! Two points
    real(KIND=wp), DIMENSION(3) :: p1, p2, p3, ser

    call eudist(p1, p2, ser(1))
    call eudist(p1, p3, ser(2))
    call eudist(p2, p2, ser(3))
end subroutine get_ser

subroutine get_delta_ser(p1, p2, p3, delta, ser)
    use kinds, ONLY: wp => dp
    implicit none
    ! points
    real(KIND=wp), DIMENSION(3) :: p1, p2, p3
    ! return n x d (3) vector 
    real(KIND=wp), DIMENSION(3, 3) :: ser
    real(KIND=wp) :: delta
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
    real(KIND=wp) :: delta, dist
    ! Dimension which to add delta
    integer :: dim
    p1(dim) = p1(dim) + delta
    p2(dim) = p2(dim) + delta
    call eudist(p1, p2, dist)
    ! Not divide by delta to get "derivative"
    ! dist = dist / delta
end subroutine eudist_with_delta

!    call compute_force( x(1, :), x(3, :), delta f_AC)
subroutine compute_force(p1, p2, p3, delta, force)
    use kinds, ONLY: wp => dp
    implicit none
    ! particles
    real(KIND=wp), DIMENSION(3) :: p1, p2, p3
    ! distances between particles
    real(KIND=wp) :: d_AB, d_AC, d_BC
    ! Vector of forces betweem each particle
    real(KIND=wp), DIMENSION(3) :: force
    !
    real(KIND=wp), DIMENSION(3) :: ser, er, der
    real(KIND=wp), DIMENSION(3) :: delta_ser, delta_er, delta_der
    ! delta to add to dimension, and distance to return
    real(KIND=wp) :: delta, dist
    ! for looping
    integer :: dim, atom

    ! get der first for actual location
    ! First we get interatomic distances
    call eudist(p1, p2, d_AB)
    call eudist(p1, p3, d_AC)
    call eudist(p2, p3, d_BC)

    ser = (/d_AB, d_AC, d_BC/)
    call jpca15(ser, er, der)

    ! now ser for location + delta (just in x direction)
    ! First we get interatomic distances of location plus delta
    call eudist_with_delta(p1, p2, delta, 1, d_AB)
    call eudist_with_delta(p1, p3, delta, 1, d_AC)
    call eudist_with_delta(p2, p3, delta, 1, d_BC)

    delta_ser = (/d_AB, d_AC, d_BC/)
    call jpca15(delta_ser, delta_er, delta_der)

    do atom = 1, 3, 1
        force(atom) = (delta_ser(atom) - ser(atom)) / delta
        ! print *, "force(atom): ", atom, force(atom)
    end do
end subroutine compute_force
