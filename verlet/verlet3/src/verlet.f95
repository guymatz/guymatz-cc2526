module verlet

    ! use kinds, ONLY: wp => dp
    ! implicit none
    ! private
    ! public :: eudist, eudist_with_delta, compute_force, get_ser, get_delta_ser

contains

    real(KIND=wp) function eudist(p1, p2) result(dist)
        use kinds, ONLY: wp => dp
        implicit none
        ! Two points
        real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2
        ! Distance to return
        !real(KIND=wp) :: dist
        dist = SQRT((p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2)
        !print *, "poop"
    end function eudist

    subroutine get_ser(p1, p2, p3, ser)
        !  Returns the distances between pairs of points
        use kinds, ONLY: wp => dp
        implicit none
        ! Three points
        real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2, p3
        real(KIND=wp), DIMENSION(3), intent(out) :: ser
        ! real(KIND=wp) :: eudist

        !ser(1) = eudist(p1, p2)
        ser(1) = eudist(p1, p2)
        ser(2) = eudist(p1, p3)
        ser(3) = eudist(p2, p3)
    end subroutine get_ser

    subroutine get_delta_ser(p1, p2, p3, delta, ser)
        !  Returns the distances between pairs of points, with a delta added in each dimension
        use kinds, ONLY: wp => dp
        implicit none
        ! points
        real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2, p3
        real(KIND=wp), intent(in) :: delta
        ! Return n x d (3) vector
        real(KIND=wp), DIMENSION(3, 3), intent(out) :: ser
        ! for looping
        integer :: dimn

        do dimn = 1, 3, 1
            ser(1, dimn) = eudist_with_delta(p1, p2, delta, dimn)
            ser(2, dimn) = eudist_with_delta(p1, p3, delta, dimn)
            ser(3, dimn) = eudist_with_delta(p2, p3, delta, dimn)
        end do
    end subroutine get_delta_ser

    real(KIND=wp) function eudist_with_delta(p1, p2, delta, dimn) result(dist)
        use kinds, ONLY: wp => dp
        implicit none
        ! Two points
        real(KIND=wp), DIMENSION(3), intent(in) :: p1, p2
        ! delta to add to dimension, and distance to return
        real(KIND=wp), intent(in) :: delta
        ! Dimension which to add delta
        integer, intent(in) :: dimn
        ! distance to return
        !real(KIND=wp), intent(out) :: dist
        !!!!!   NOT Adding delta to second atom.  Only the first
        dist = SQRT((p1(dimn) + delta - p2(dimn))**2)
        ! print *, "Dim:", dimn, "p1:", p1(dimn), "p2:", p2(dimn), " -> ", dist
    end function eudist_with_delta

    ! points = x,y,z coords of each point
    subroutine compute_force(points, delta, forces)
        use kinds, ONLY: wp => dp
        implicit none
        ! IN
        real(KIND=wp), DIMENSION(3, 3), intent(in) :: points
        real(KIND=wp), intent(in) :: delta
        ! OUT
        real(KIND=wp), DIMENSION(3, 3), intent(out) :: forces
        ! We will compute `ser` for input to jpca15, and stores the results in der
        real(KIND=wp), DIMENSION(3) :: ser, der
        ! jpca15 also returns er
        real(KIND=wp) :: er
        ! single_force is used to store the force on an atom in one direction/dimension
        !! real(KIND=wp) :: single_force
        ! We will compute `delta_ser` for input to jpca15, and stores the results in delta_der
        real(KIND=wp), DIMENSION(3) :: delta_ser, delta_der
        ! jpca15 also returns delta_er
        real(KIND=wp) :: delta_er
        ! d_ is used to store the euclidean distances between atoms
        real(KIND=wp), DIMENSION(3) :: d_
        ! delta_d_ is used to store the euclidean distances + delta between atoms
        ! It's 2D so that we can store distances between atoms in all dimensions
        real(KIND=wp), DIMENSION(3, 3) :: delta_d_
        ! for looping
        integer :: dimn, atom_i

        ! Step 1 of calculating forces: - get euclidean distances between points
        d_(1) = eudist(points(1, :), points(2, :))
        d_(2) = eudist(points(1, :), points(3, :))
        d_(3) = eudist(points(2, :), points(3, :))
        ! Step 2 of calculating forces: - pass distances as to jpca15 as `ser`
        ser = (/d_(1), d_(2), d_(3)/)
        call jpca15(ser, er, der)
        ! Step 2 of calculating forces: - and get back er & der.  We dont use `er`

        ! Step 3 of calculating forces: - get euclidean distances between points + delta
        do dimn = 1, 3, 1
            delta_d_(1, dimn) = eudist_with_delta(points(1, :), points(2, :), delta, dimn)
            delta_d_(2, dimn) = eudist_with_delta(points(1, :), points(3, :), delta, dimn)
            delta_d_(3, dimn) = eudist_with_delta(points(2, :), points(3, :), delta, dimn)
        end do

        ! Step 3.5: I don't think I need to do this, but let's initialize the forces to 0
        do dimn = 1, 3, 1
            do atom_i = 1, 3, 1
                forces(atom_i, dimn) = 0.0_wp
            end do
        end do

        ! Step 4: Update the nx3 `forces` array with the force on each atom in each dimension
        do atom_i = 1, 3, 1
            do dimn = 1, 3, 1
                ! print *, ""
                ! print *, "Delta:", delta
                delta_ser = (/delta_d_(1, dimn), delta_d_(2, dimn), delta_d_(3, dimn)/)
                call jpca15(delta_ser, delta_er, delta_der)
                ! single_force = (delta_der(dimn) - der(dimn)) / delta
                ! forces(atom_i, dimn) = single_force
                forces(atom_i, dimn) = (delta_der(dimn) - der(dimn)) / delta

                ! print *, "SER:", ser
                ! print *, "DER:", der
                ! print *, "delta_SER:", delta_ser
                ! print *, "delta_DER:", delta_der
                ! print *, "Atom #:", atom_i, ", Dimension:", dimn, ", Force:", forces(atom_i, dimn) 

            end do
        end do
        ! DONE!
    end subroutine compute_force

end module verlet
