program main
! 1verbose set fdm?

! Var definitions, etc
    ! For printing to STDERR
    use,intrinsic :: iso_fortran_env, only : stderr=>ERROR_UNIT
    use kinds, ONLY: wp => dp
    use verlet
    implicit none

! vars
    ! for looping by atoms & dimension
    integer :: atom_num, dimn

    LOGICAL :: OK = .FALSE.
    ! For requesting xyz file output
    LOGICAL :: XYZ = .FALSE.
    CHARACTER(len=32) :: arg
    CHARACTER(len=32) :: file_name = "atoms.dat"
    !integer :: arg_len
    !integer :: status
    real(KIND=wp), DIMENSION(:, :), ALLOCATABLE :: x, v, f, fnext
    real(KIND=wp), DIMENSION(:), ALLOCATABLE :: mass
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
    integer :: nk, nk_cli=0
    !real :: sigma, epsilon
    real(KIND=wp) :: tau, tau_cli=0.0_wp
    !real :: tmp
    ! interatomic distances
    !real(KIND=wp) :: d_AB, d_AC, d_BC
    ! Delta for interatomic distances
    !real(KIND=wp), DIMENSION(3) :: dd_AB, dd_AC, dd_BC
    real(KIND=wp) :: delta = 0.1, delta_cli = 0
    ! For computing jPCA with delta
    !real(KIND=wp), DIMENSION(3) :: dd_ser, dd_er, dd_der
    !integer :: num_rows
    integer :: num_atoms
    ! locations and velocites for atom
    real(KIND=wp) :: ax, ay, az, vx, vy, vz
    ! for storing intermediate values & looping
    integer :: i, k
    !real (KIND = wp), DIMENSION(7) :: p_a, p_b ! our two particles

! Process command-line args
    i = 0
    DO
        i = i + 1
        IF (i .gt. command_argument_count()) exit

        CALL get_command_argument(i, arg)
        ! Delta - defaults to 0.01 (see above)
        IF (arg == "-d") THEN
            ! delta
            i = i + 1
            CALL get_command_argument(i, arg)
            write(stderr,*) "arg -d:", arg
            read (arg, '(f33.32)') delta_cli
            IF (delta_cli == 0.0) THEN
                print *, "Delta too small!"
                print '(A,I0)', "main +", __LINE__
                STOP __LINE__ - 1
            END IF
        ELSE IF (arg == "-s") THEN
            ! Num steps - defaults to 0 (see above): Should be in data file
            i = i + 1
            CALL get_command_argument(i, arg)
            read (arg, '(I5)') nk_cli
        ELSE IF (arg == "-t") THEN
            ! tau
            i = i + 1
            CALL get_command_argument(i, arg)
            read (arg, '(f1.2)') tau_cli
        ELSE IF (arg == "-f") THEN
            ! filename
            i = i + 1
            CALL get_command_argument(i, arg)
            file_name = trim(arg)
            INQUIRE (FILE=file_name, EXIST=OK)
            if (.NOT. OK) THEN
                print *, "ERROR!!  File does not exist: ", file_name
                print '(A,I0)', "main +", __LINE__
                STOP __LINE__ - 1
            END IF
        ELSE IF (arg == "-x") THEN
            ! Print xyz format
            XYZ = .TRUE.
        ELSE
            print *, "Arg used: ", arg
            print *, "Usage: verlet3 [ -h ] [ -f data_file ]  [ -d delta ] [ -s steps ] [-x]"
            print *, "DEFAULTS:"
            print '(A, A)', "    data_file - ", file_name
            print '(A, F0.9)', "        delta - ", delta
            print '(A, F0.9)', "        tau - ", tau
            print '(A, I0)', "        steps - ", nk
            print '(L1)', "     xyz file - ", XYZ
            print '(A,I0)', "main +", __LINE__
            STOP __LINE__ - 1
        END IF
    END DO

! Read atoms, etc from data file
    open (UNIT=11, FILE=file_name, STATUS="old", ACTION="read")
    read (unit=11, FMT=*) nk, tau
    read (unit=11, FMT=*) delta
    read (unit=11, FMT=*) num_atoms
    !print *, nk, tau, sigma, epsilon, num_atoms
    write(stderr,*) "Number of atoms:", num_atoms

    ! Let the command-line `steps` override nk, if it's set
    if (nk_cli > 0) then
        nk = nk_cli
    end if
    ! same for delta
    if (delta_cli > 0) then
        delta = delta_cli
    end if
    ! nd tau
    if (tau_cli > 0) then
        tau = tau_cli
    end if

    write(stderr,*) "Will use delta: ", delta
    write(stderr,*) "Will use file: ", file_name
    write(stderr,*) "Will use num steps: ", nk
    write(stderr,*) "Will use tau: ", tau
    write(stderr,*) "Will print XYZ file: ", XYZ

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
    allocate (mass(num_atoms))

    ! read in info for particles
    do i = 1, num_atoms, 1
        read (unit=11, FMT=*) mass(i), ax, ay, az, vx, vy, vz
        x(i, :) = (/ax, ay, az/)
        v(i, :) = (/vx, vy, vz/)
        write(stderr,*) 'Particle', i, ': '
        write(stderr,*) '  Mass:', mass(i)
        write(stderr,*) '  Starting Position:', x(i, :)
        write(stderr,*) '  Initial Velocity: ', v(i, :)
    end do
    close (unit=11)

! print initial output for XYZ data file -
    ! https://en.wikipedia.org/wiki/XYZ_file_format
    if (XYZ) then
        print *, num_atoms
        print *, "Initial Positions"
        print *, "atom1", x(1, :)
        print *, "atom2", x(2, :)
        print *, "atom3", x(3, :)
    end if
! Initial Force for particles
    call compute_force(x, delta, f)
    write(stderr,*) " ***** Initial Positions / Forces: "
    do i = 1, 3, 1
        write(stderr,*) "  Atom:", i
        write(stderr,*) "     x: ", x(i, :)
        write(stderr,*) "     f: ", f(i, :)
    end do

  ! Iterate!
    ! write(stderr,*) "  nk:", nk
    do k = 1, nk, 1
        ! Calculate x^{(a)}_{k+1}
        ! Calculate new position

        do atom_num = 1, 3, 1
            print *, " ***** x before:", atom_num, x(atom_num, :) 
        end do
        do atom_num = 1, 3, 1
            do dimn = 1, 3, 1

                ! print *, "NUMS: ", atom_num, dimn
                ! print *, "      x: ", x(atom_num, dimn)
                ! print *, "      t: ", tau
                ! print *, "      v: ", v(atom_num, dimn)
                ! print *, "      f: ", f(atom_num, dimn)
                ! print *, "      m: ", mass(atom_num)
                ! print *, "     l1: ", tau * v(atom_num, dimn)
                ! print *, "     l2: ", x(atom_num, dimn) + tau * v(atom_num, dimn)
                ! print *, "   l3-1: ", f(atom_num, dimn)
                ! print *, "   l3-2: ", 2 * mass(atom_num)
                ! print *, "     l3: ", f(atom_num, dimn) / (2 * mass(atom_num)) * tau**2

                x(atom_num, dimn) = x(atom_num, dimn) + tau * v(atom_num, dimn) + &
                                    (f(atom_num, dimn) / (2 * mass(atom_num))) * tau**2
                print *, "  new x: ", x(atom_num, dimn)
            end do
        end do
        do atom_num = 1, 3, 1
            print *, " ***** x AFTER:", atom_num, x(atom_num, :) 
        end do
        ! calculate fnext
        call compute_force(x, delta, fnext)
        do atom_num = 1, 3, 1
            print *, "POOP fn:", atom_num, fnext(atom_num, :) 
        end do
        ! Calculate new velocity
        ! A
        do atom_num = 1, 3, 1
            do dimn = 1, 3, 1
                v(atom_num, dimn) = v(atom_num, dimn) + &
                                    tau / (2 * mass(atom_num)) * (f(atom_num, dimn) + &
                                    fnext(atom_num, dimn))
            end do
        end do

        if (XYZ) then
            print *, num_atoms
            print *, "step:", k
            print *, "atom1", x(1, :)
            print *, "atom2", x(2, :)
            print *, "atom3", x(3, :)
        end if

        ! Assign f = fnext
        do atom_num = 1, 3, 1
          do dimn = 1, 3, 1
            f(atom_num, dimn) = fnext(atom_num, dimn)
            ! print *, atom_num, dimn, f(atom_num, dimn) 
          end do
        end do

    end do

! Cleanup
    deallocate (x)
    deallocate (v)
    deallocate (f)
    deallocate (fnext)
    deallocate (mass)

end program main
