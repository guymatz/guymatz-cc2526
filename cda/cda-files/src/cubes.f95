MODULE cubes
  USE kinds, ONLY: wp => dp
  IMPLICIT NONE

  TYPE, PUBLIC :: cube
    PRIVATE
    CHARACTER (LEN=72) :: data_title = 'Untitled'
    CHARACTER (LEN=72) :: data_desc = 'No Description'
    real (KIND=wp) :: xmin, ymin, zmin, dx, dy, dz
    integer :: nx, ny, nz, natom = 0
    integer, DIMENSION(:), POINTER :: zahl
    real (KIND=wp), DIMENSION(:), POINTER :: chrg, x, y, z
    real (KIND=wp), DIMENSION(:), POINTER :: array
    ! fill
  END TYPE cube

  PUBLIC :: cube_get, &
            cube_add, &
            cube_sub, &
            cube_int, &
            cube_cdz, &
            cube_del, &
            cube_gen_charge_redistribution, &
            cube_print

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE cube_add
  END INTERFACE
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE cube_sub
  END INTERFACE


  CONTAINS


  SUBROUTINE cube_get (mycube, infile)
    CHARACTER(LEN=*), INTENT(IN) :: infile
    TYPE (cube), INTENT(OUT) :: mycube
    ! access attributes as follows:
    ! mycube%data_title = ...
    ! temp var for unwannted data
    real :: tmp
    integer :: i, x, y, z
    OPEN (UNIT=11, FILE=infile, STATUS="old", ACTION="read")
    read(11, '(A)') mycube%data_title
    read(11, '(A)') mycube%data_desc
    read(11, FMT=*) mycube%natom, mycube%xmin, mycube%ymin, mycube%zmin
    read(11, FMT=*) mycube%nx, mycube%dx, tmp, tmp
    read(11, FMT=*) mycube%ny, tmp, mycube%dy, tmp
    read(11, FMT=*) mycube%nz, tmp, tmp, mycube%dz
    ! 2 x 5 array - data looks like
    ! ZNatoms, chargeNatoms, xNatoms, yNatoms, zNatoms
    allocate(mycube%zahl(mycube%natom))
    allocate(mycube%chrg(mycube%natom))
    allocate(mycube%x(mycube%natom))
    allocate(mycube%y(mycube%natom))
    allocate(mycube%z(mycube%natom))
    do i = 1, mycube%natom, 1
      read(11, FMT=*) mycube%zahl(i), mycube%chrg(i), mycube%x(i), mycube%y(i), mycube%z(i)
    end do
    ! slurp the rest of the file into c%array
    allocate(mycube%array(mycube%nx * mycube%ny * mycube%nz))
    read(11, FMT=*) mycube%array
  END SUBROUTINE cube_get

  subroutine cube_print(c)
    TYPE(cube), INTENT(IN) :: c
    integer :: i, ctr = 4
    print *, 'Title: ', c%data_title
    print *, 'Desc.: ', c%data_desc
    print *, 'Num Atoms: ', c%natom
    print *, 'Atoms:'
    do i = 1, c%natom, 1
        print *, c%zahl(i), c%chrg(i), c%x(i), c%y(i), c%z(i)
    end do
    print '("First ", I0, " of Array:")', ctr 
    do i = 1, ctr, 1
        print *, c%array(i)
    end do
  end subroutine cube_print

  FUNCTION cube_add (mycube1, mycube2)
    TYPE(cube) :: cube_add
    TYPE(cube), INTENT(IN) :: mycube1, mycube2
    allocate(cube_add%array(size(mycube1%array)))
    cube_add%array = mycube1%array + mycube2%array
  END FUNCTION cube_add


  FUNCTION cube_sub (mycube1, mycube2)
    TYPE(cube) :: cube_sub
    TYPE(cube), INTENT(IN) :: mycube1, mycube2
    allocate(cube_sub%array(size(mycube1%array)))
    cube_sub%array = mycube1%array - mycube2%array
  END FUNCTION cube_sub


  FUNCTION cube_int (mycube)
    REAL (KIND=wp) :: cube_int
    TYPE (cube), INTENT(IN) :: mycube
    ! ...
  END FUNCTION cube_int


  SUBROUTINE cube_del (mycube)
    TYPE (cube), INTENT(IN) :: mycube
    ! ...
  END SUBROUTINE cube_del

  ! Take in input a cube variable and returning a one-dimensional array
  ! containing the values of Δq
  SUBROUTINE cube_cdz (mycube)
    TYPE (cube), INTENT(IN) :: mycube
    ! ...
  END SUBROUTINE cube_cdz

  ! Take 4 cubes - a reference, A, B, & AB  - and set the reference%array
  ! to AB%array - A%array - B%array
  SUBROUTINE cube_gen_charge_redistribution (delta, AB, A, B)
    TYPE (cube), INTENT(IN) :: AB, A, B
    TYPE (cube), INTENT(OUT) :: delta
    allocate(delta%array(size(AB%array)))
    delta%array = AB%array - A%array - B%array
  END SUBROUTINE cube_gen_charge_redistribution

END MODULE cubes
