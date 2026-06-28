program test
! preamble
    use kinds, ONLY: wp => dp
    use verlet
    implicit none

! variables
    real(KIND=wp), DIMENSION(:, :), ALLOCATABLE :: x, v, f, fnext, points, forces
    real(KIND=wp), DIMENSION(:), ALLOCATABLE :: mass
    real(KIND=wp), DIMENSION(3) :: distance_with_delta
    integer :: nk ! num interations
    integer :: n = 3  ! num atoms
    real(KIND=wp) :: tau, delta
    real(KIND=wp) :: dist
    allocate (x(n, 3))
    allocate (v(n, 3))
    allocate (f(n, 3))
    allocate (fnext(n, 3))
    allocate (mass(n))
    allocate (points(n, 3))
    allocate (forces(n, 3))

! Initialization
    nk = 6000                                     ! nk
    tau = 0.2                                     ! tau
    delta = 0.1                                          ! delta
    ! atom A
    mass(1) = 1.0080
    x(1, :) = (/-10.0, 0.0, 0.0/)          ! x, y, z
    v(1, :) = (/1.0, 0.0, 0.0/)        ! vx, vy, vz
    ! atom B
    mass(2) = 1.0080                               ! m, x, y, z, vx, vy, vz
    x(2, :) = (/10.0, 0.0, 0.0/)           ! m, x, y, z, vx, vy, vz
    v(2, :) = (/0.0, 0.0, 0.0/)        ! m, x, y, z, vx, vy, vz
    ! atom C
    mass(3) = 1.0080                               ! m, x, y, z, vx, vy, vz
    x(3, :) = (/10.74, 0.0, 0.0/)           ! m, x, y, z, vx, vy, vz
    v(3, :) = (/0.0, 0.0, 0.0/)        ! m, x, y, z, vx, vy, vz

! meat
    print *, "LOCATION ==========="
    print *, "Location atom A: ", x(1, :)
    print *, "Location atom B: ", x(2, :)
    print *, "Location atom C: ", x(3, :)
    print *, " TEST eudist distance ====================="
    print *, 'eudist: A - B', eudist(x(1, :), x(2, :))
    print *, 'eudist: A - C', eudist(x(1, :), x(3, :))
    print *, 'eudist: B - C', eudist(x(2, :), x(3, :))
    print *, " TEST distance w/ delta ========      x                        y              &
&          z"
    distance_with_delta = eudist_with_delta(x(1, :), x(2, :), delta, 1)
    print *, 'dist w/ delta: A - B, d x', dist
    distance_with_delta = eudist_with_delta(x(1, :), x(2, :), delta, 2)
    print *, 'dist w/ delta: A - B, d y', distance_with_delta
    distance_with_delta = eudist_with_delta(x(1, :), x(2, :), delta, 3)
    print *, 'dist w/ delta: A - B, d z', distance_with_delta
    print *, ""
    distance_with_delta = eudist_with_delta(x(1, :), x(3, :), delta, 1)
    print *, 'dist w/ delta: A - C, d x', distance_with_delta
    distance_with_delta = eudist_with_delta(x(1, :), x(3, :), delta, 2)
    print *, 'dist w/ delta: A - C, d y', distance_with_delta
    distance_with_delta = eudist_with_delta(x(1, :), x(3, :), delta, 3)
    print *, 'dist w/ delta: A - C, d z', distance_with_delta
    print *, ""
    distance_with_delta = eudist_with_delta(x(2, :), x(3, :), delta, 1)
    print *, 'dist w/ delta: B - C, d x', distance_with_delta
    distance_with_delta = eudist_with_delta(x(2, :), x(3, :), delta, 2)
    print *, 'dist w/ delta: B - C, d y', distance_with_delta
    distance_with_delta = eudist_with_delta(x(2, :), x(3, :), delta, 3)
    print *, 'dist w/ delta: B - C, d z', distance_with_delta
    print *, ""

    print *, " TEST compute_force ========      x                        y                     z"
    PRINT *, forces
    call compute_force(x, delta, forces)
    print *, " Forces ========      x                        y                     z"
    print *, 'forces: A - B - C', forces

! clean up
    deallocate (x)
    deallocate (v)
    deallocate (f)
    deallocate (fnext)
    deallocate (mass)

end program
