! subroutine eudist(p1, p2, dist)
!     use kinds, ONLY: wp => dp
!     implicit none
!     ! Two points
!     real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2
!     ! Distance to return
!     real(KIND=wp), intent(out) :: dist
!     dist = SQRT((p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2)
! end subroutine eudist

pure function eudist(p1, p2) result(dist)
    use kinds, ONLY: wp => dp
    implicit none
    ! Two points
    real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2
    ! Distance to return
    real(KIND=wp) :: dist
    dist = SQRT((p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2)
end function eudist

subroutine get_ser(p1, p2, p3, ser)
    use kinds, ONLY: wp => dp
    implicit none
    ! Three points
    real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2, p3
    real(KIND=wp), DIMENSION(3), intent(out) :: ser
    real(KIND=wp) :: eudist

    !ser(1) = eudist(p1, p2)
    ser(1) = eudist(p1, p2)
    ser(2) = eudist(p1, p3)
    ser(3) = eudist(p2, p2)
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
    !!!!!   NOT Adding delta to second atom.  Only the first
    dist = SQRT((p1(dim) + delta - p2(dim))**2)
end subroutine eudist_with_delta

!    call compute_force(x, delta f_AC)
! points = x,y,z coords of each point
subroutine compute_force(points, delta, forces)
    use kinds, ONLY: wp => dp
    implicit none
    !integer :: n
    ! particles
    real(KIND=wp), DIMENSION(3, 3), intent(in) :: points
    real(KIND=wp), intent(in) :: delta
    ! distances between particles
    !real(KIND=wp) :: d_AB, d_AC, d_BC
    ! Vector of forces betweem each particle
    real(KIND=wp), DIMENSION(3, 3), intent(inout) :: forces
    !
    real(KIND=wp), DIMENSION(3) :: ser, der
    real(KIND=wp) :: er, single_force
    real(KIND=wp), DIMENSION(3) :: delta_ser, delta_der
    real(KIND=wp) :: delta_er
    ! n*(n-1)/2 gets number of pairs of particles
    real(KIND=wp), DIMENSION(3) :: d_
    ! n*(n-1)/2 gets number of pairs of particles in each direction
    real(KIND=wp), DIMENSION(3, 3) :: delta_d_
    ! delta to add to dimension, and distance to return
    ! for looping
    !integer :: dim, atom_i, atom_j
    integer :: dim, atom_i
    real(KIND=wp) :: eudist

    print *, "p1: ", points(1, :)
    ! get der first for actual location
    ! First we get interatomic distances
    d_(1) = eudist(points(1, :), points(2, :))
    d_(2) = eudist(points(1, :), points(3, :))
    d_(3) = eudist(points(2, :), points(3, :))
    ! now ser for location + delta

    ser = (/d_(1), d_(2), d_(3)/)
    call jpca15(ser, er, der)
    print *, "p2: ", points(2, :)

    do dim = 1, 3, 1
        call eudist_with_delta(points(1, :), points(2, :), delta, dim, delta_d_(1, dim))
        call eudist_with_delta(points(1, :), points(3, :), delta, dim, delta_d_(2, dim))
        call eudist_with_delta(points(2, :), points(3, :), delta, dim, delta_d_(3, dim))
    end do

    PRINT *, forces
    do dim = 1, 3, 1
        do atom_i = 1, 3, 1
            PRINT *, forces(atom_i, dim), atom_i, dim
            forces(atom_i, dim) = 0.0_wp
        end do
    end do
    STOP 124

    print *, "p3: ", points(3, :)
    do dim = 1, 3, 1
        do atom_i = 1, 3, 1
            print *, "dim / atom: ", dim, atom_i
            delta_ser = (/delta_d_(1, dim), delta_d_(2, dim), delta_d_(3, dim)/)
            print *, "delta_ser: ", delta_ser
            call jpca15(delta_ser, delta_er, delta_der)
            print *, "delta_der: ", delta_der
            print *, "der: ", der
            print *, "delta: ", delta
            print *, "derivative:", (delta_der(dim) - der(dim)) / delta
            print *, "size of force: ", size(forces)
            print *, "size of force(1): ", size(forces(1, :))
            print *, "size of force(2): ", size(forces(2, :))
            print *, "size of force(3): ", size(forces(3, :))
            single_force = (delta_der(dim) - der(dim)) / delta
            print *, "Force: ", single_force
            print *, "Force 1, 1: ", forces(atom_i, dim)
            forces(atom_i, dim) = single_force
        end do
    end do

    print *, "d_(1): ", d_(1)
    print *, "d_(2): ", d_(2)
    print *, "d_(3): ", d_(3)

    print *, "delta_d_(1): ", delta_d_(1, :)
    print *, "delta_d_(2): ", delta_d_(2, :)
    print *, "delta_d_(3): ", delta_d_(3, :)

    print *, "force_(1): ", forces(1, :)
    print *, "force_(2): ", forces(2, :)
    print *, "force_(3): ", forces(3, :)
    stop 3

    ! ser = (/d_AB, d_AC, d_BC/)
    ! call jpca15(ser, er, der)
    ! call eudist_with_delta(p1, p3, delta, 1, d_AC)
    ! call eudist_with_delta(p2, p3, delta, 1, d_BC)
    ! delta_ser = (/d_AB, d_AC, d_BC/)
    ! call jpca15(delta_ser, delta_er, delta_der)
end subroutine compute_force
