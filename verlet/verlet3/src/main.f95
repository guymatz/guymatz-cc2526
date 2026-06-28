program main

! Var definitions, etc
    use kinds, ONLY: wp => dp
    use verlet
    implicit none

! vars
    ! for looping by atoms & dimension
    integer :: atom_num, dim

    LOGICAL :: OK
    CHARACTER(len=32) :: arg
    CHARACTER(len=32) :: file_name = "atoms.dat"
    !integer :: arg_len
    !integer :: status
    real(KIND=wp), DIMENSION(:, :), ALLOCATABLE :: x, v, f, fnext, mass
    ! ser, er & der are for parameters to jpca15 function
    ! INPUT
    !   ser: a vector with the three interatomic distances (AB, AC, and BC)
    ! OUTPUT
    !   er: potential energy (in eV)
    !   der: vector of the derivatives of the potential with
    !        respect to the interatomic distances (AB, AC, and BC) (in bohr)
    !real(KIND=wp), DIMENSION(3) :: d_AB_ser, d_AB_er, d_AB_der
    !real(KIND=wp), DIMENSION(3) :: d_AC_ser, d_AC_er, d_AC_der
    !real(KIND=wp), DIMENSION(3) :: d_BC_ser, d_BC_er, d_BC_der
    ! Force vectors (in the x-direction only)
    !real(KIND=wp), DIMENSION(3) :: f_AB, f_AC, f_BC, force
    ! for reading in atomic data from file
    integer :: nk
    !real :: sigma, epsilon
    real(KIND=wp) :: tau
    !real :: tmp
    ! interatomic distances
    !real(KIND=wp) :: d_AB, d_AC, d_BC
    ! Delta for interatomic distances
    !real(KIND=wp), DIMENSION(3) :: dd_AB, dd_AC, dd_BC
    real(KIND=wp) :: delta = 0.1, cli_delta = 0
    ! For computing jPCA with delta
    !real(KIND=wp), DIMENSION(3) :: dd_ser, dd_er, dd_der
    !integer :: num_rows
    integer :: num_atoms
    ! locations and velocites for atom
    real(KIND=wp) :: ax, ay, az, vx, vy, vz
    ! for storing intermediate values & looping
    integer :: i, k
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
                print *, "main +67"
                STOP 67
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
                print *, "main +108"
                STOP 108
            END IF
        ELSE
            print *, "Usage: verlet3 [ -h ] [ -f data_file ]  [ -d delta ] [ -s steps ]"
            print *, "DEFAULTS:"
            print '(A, A)', "    data_file - ", file_name
            print '(A, F0.9)', "        delta - ", delta
            print '(A, I0)', "        steps - ", steps
            print *, "main +117"
            STOP 117
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
    call compute_force(x, delta, f)

    ! Iterate!
    do k = 1, nk, 1
        ! Calculate x^{(a)}_{k+1}
        do atom_num = 1, num_atoms, 1
            do dim = 1, 3, 1 ! dimensions
                ! Caclculate new position
                x(atom_num, dim) = x(atom_num, dim) + tau * v(atom_num, dim) &
                                   + tau**2 * f(atom_num, dim) / (2 * mass(atom_num, dim))
                ! calculate fnext
                call compute_force(x, delta, f)
                ! Calculate new velocity
                v(atom_num, dim) = v(atom_num, dim) + tau / (2 * mass(atom_num, dim)) * &
                                   (f(atom_num, dim) + fnext(atom_num, dim))
                ! Assign f = fnext
                f(atom_num, dim) = fnext(atom_num, dim)
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
    deallocate (mass)

end program main
