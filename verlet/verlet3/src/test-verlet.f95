program test
! preamble
    use kinds, ONLY: wp => dp
    use verlet
    implicit none

! variables
    real(KIND=wp), DIMENSION(:, :), ALLOCATABLE :: x, v, f, fnext, points, forces
    real(KIND=wp), DIMENSION(:), ALLOCATABLE :: mass
    real(KIND=wp), DIMENSION(:), ALLOCATABLE :: ser
    real(KIND=wp), DIMENSION(:, :), ALLOCATABLE :: ser_delta
    !real(KIND=wp), DIMENSION(3) :: distance_with_delta
    real(KIND=wp) :: distance_with_delta
    integer :: nk ! num interations
    integer :: n = 3  ! num atoms
    integer :: i  ! for loops!
    real(KIND=wp) :: tau, delta
    !real(KIND=wp) :: dist
    allocate (x(n, 3))
    allocate (v(n, 3))
    allocate (f(n, 3))
    allocate (fnext(n, 3))
    allocate (mass(n))
    allocate (points(n, 3))
    allocate (ser(n))
    allocate (forces(n, 3))
    allocate (ser_delta(n, 3))

! Initialization
    nk = 6000                                     ! nk
    tau = 0.2                                     ! tau
    delta = 0.5                                          ! delta
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

!! TESTS

! Initial locations & Distances
    print *, "=========== Initial LOCATION"
    print *, "Atom A: ", x(1, :)
    print *, "Atom B: ", x(2, :)
    print *, "Atom C: ", x(3, :)
    print *, ""
    print *, "=========== TEST eudist: Gets Euclidean distance between two atoms"
    print *, 'eudist: A - B: ', eudist(x(1, :), x(2, :))
    print *, 'eudist: A - C: ', eudist(x(1, :), x(3, :))
    print *, 'eudist: B - C: ', eudist(x(2, :), x(3, :))
    print *, ""
    print *, "===========  TEST eudist_with_delta: Gets distance between two atoms after adding a small delta in each dimension" 
    print *, "delta = ", delta
    print *, "                        p1                            p2                       d"

! Initial Distances w/ Delta
    distance_with_delta = eudist_with_delta(x(1, :), x(2, :), delta, 1)
    print *, 'dist w/ delta: A', x(1, 1), ' - B:', x(2, 1), ' -> d x', distance_with_delta
    distance_with_delta = eudist_with_delta(x(1, :), x(2, :), delta, 2)
    print *, 'dist w/ delta: A', x(1, 2), ' - B:', x(2, 2), ' -> d y', distance_with_delta
    distance_with_delta = eudist_with_delta(x(1, :), x(2, :), delta, 3)
    print *, 'dist w/ delta: A', x(1, 3), ' - B:', x(2, 3), ' -> d z', distance_with_delta
    print *, ""
    distance_with_delta = eudist_with_delta(x(1, :), x(3, :), delta, 1)
    print *, 'dist w/ delta: A', x(1, 1), ' - C:', x(3, 1), ' -> d x', distance_with_delta
    distance_with_delta = eudist_with_delta(x(1, :), x(3, :), delta, 2)
    print *, 'dist w/ delta: A', x(1, 2), ' - C:', x(3, 2), ' -> d y', distance_with_delta
    distance_with_delta = eudist_with_delta(x(1, :), x(3, :), delta, 3)
    print *, 'dist w/ delta: A', x(1, 3), ' - C:', x(3, 3), ' -> d z', distance_with_delta
    print *, ""
    distance_with_delta = eudist_with_delta(x(2, :), x(3, :), delta, 1)
    print *, 'dist w/ delta: B', x(2, 1), ' - C:', x(3, 1), ' -> d x', distance_with_delta
    distance_with_delta = eudist_with_delta(x(2, :), x(3, :), delta, 2)
    print *, 'dist w/ delta: B', x(2, 2), ' - C:', x(3, 2), ' -> d y', distance_with_delta
    distance_with_delta = eudist_with_delta(x(2, :), x(3, :), delta, 3)
    print *, 'dist w/ delta: B', x(2, 3), ' - C:', x(3, 3), ' -> d z', distance_with_delta

! Test get_ser
    ! subroutine get_ser(p1, p2, p3, ser)
    print *, ""
    print *, "======== Test get_ser: Gets distances between 3 atoms as (AB, AC, BC)"
    CALL get_ser(x(1, :), x(2, :), x(3, :), ser)
    do i = 1, 3, 1
        print *, "Distance :", i, ":", ser(i)
    end do

! Test get_delta_ser
    ! subroutine get_delta_ser(p1, p2, p3, delta, ser)
    print *, ""
    print *, "======== Test get_delta_ser: Gets distances between 3 atoms as (AB, AC, BC) after adding small delta in each dimension"
    print *, "delta = ", delta
    print *, "                             x                             y                   z"
    CALL get_delta_ser(x(1, :), x(2, :), x(3, :), delta, ser_delta)
    do i = 1, 3, 1
        print *, "Distance :", i, ":", ser_delta(i, :)
    end do

! Test compute_force
    print *, "", "======== TEST compute_force: "
    print *, " Initial forces:  Just to show that we are starting from (0, 0, 0)"
    print *, "                             x                             y                   z"
    do i = 1, 3, 1
        PRINT *, "Atom #", i, ":", forces(i, :)
    end do

    call compute_force(x, delta, forces)
    print *, "", "======== Computed Forces: Uses jpca15 to get `der` and `delta_der` to compute  (delta_der - der) / delta"
    print *, "delta = ", delta
    print *, "                             x                             y                   z"
    do i = 1, 3, 1
        PRINT *, "Atom #", i, ":", forces(i, :)
    end do

! clean up
    deallocate (x)
    deallocate (v)
    deallocate (f)
    deallocate (fnext)
    deallocate (mass)

end program
