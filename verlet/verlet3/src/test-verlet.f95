program test
! variables
    use kinds, ONLY: wp => dp
    implicit none
    real(KIND=wp), DIMENSION(:, :), ALLOCATABLE :: x, v, f, fnext
    real(KIND=wp), DIMENSION(:), ALLOCATABLE :: mass
    real(KIND=wp) :: eudistf
    integer :: nk ! num interations
    integer :: n = 3  ! num atoms
    real(KIND=wp) :: tau, delta
    allocate (x(n, 3))
    allocate (v(n, 3))
    allocate (f(n, 3))
    allocate (fnext(n, 3))
    allocate (mass(n))

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
    print *, 'x(1, :) = ', x(1, :)
    print *, 'x(2, :) = ', x(2, :)
    print *, 'x(3, :) = ', x(3, :)

    print *, 'fnext: -2.209358416015258e-05 = '!, call jpca15((1, 3)
    print *, 'fnext: 2.209358416015258e-05 = '!, lj_out(2, 3)
    print *, 'eudistf: ', eudistf(x(1, :), x(2, :))

! clean up
    deallocate (x)
    deallocate (v)
    deallocate (f)
    deallocate (fnext)
    deallocate (mass)

end program

! real(KIND=wp) function eudistf(p1, p2)
!     use kinds, ONLY: wp => dp
!     implicit none
!     real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2
!     ! print *, "in func: eudistf before . . .  "
!     eudistf = 0.3
!     ! print *, "in func: eudistf: ", eudistf
! end function
