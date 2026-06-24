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

    !call eudist(p1, p2, ser(1))
    ser(1) = eudistf(p1, p2)
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
