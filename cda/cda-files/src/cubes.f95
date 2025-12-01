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
    do i = 1, c%natom, 1
        print *, c%zahl(i), c%chrg(i), c%x(i), c%y(i), c%z(i)
    end do
    print '("First ", I0, " of Array:")', ctr 
    do i = 1, ctr, 1
        print *, c%array(i)
    end do
  end subroutine cube_print

  FUNCTION cube_add (c1, c2)
    integer :: i
    TYPE(cube) :: cube_add
    TYPE(cube), INTENT(IN) :: c1, c2
    allocate(cube_add%array(size(c1%array)))
    allocate(cube_add%zahl(size(c1%zahl) + size(c2%zahl)))
    allocate(cube_add%chrg(size(c1%chrg) + size(c2%chrg)))
    allocate(cube_add%x(size(c1%x) + size(c2%x)))
    allocate(cube_add%y(size(c1%y) + size(c2%y)))
    allocate(cube_add%z(size(c1%z) + size(c2%z)))
    cube_add%array = c1%array + c2%array
    cube_add%data_title = c1%data_title // ' + ' // c2%data_title 
    cube_add%data_desc = c1%data_desc // ' + ' // c2%data_desc 
    !real (KIND=wp) :: xmin, ymin, zmin, dx, dy, dz
    allocate(cube_add%zahl(size(c1%zahl) + size(c2%zahl)))
    print *, 'Length of new zahl: ', size(cube_add%zahl)
    do i = 1, size(c1%zahl), 1
        cube_add%zahl(i) = c1%zahl(i)
    end do
    do i = 1, size(c2%zahl), 1
        cube_add%zahl(i + size(c1%zahl)) = c2%zahl(i)
    end do
    cube_add%chrg(i) = c1%chrg(i) + c2%chrg(i)
    cube_add%x(i) = c1%x(i) + c2%x(i)
    cube_add%y(i) = c1%y(i) + c2%y(i)
    cube_add%z(i) = c1%z(i) + c2%z(i)

    cube_add%nx = c1%nx
    cube_add%ny = c1%ny
    cube_add%nz = c1%nz
    cube_add%natom = c1%natom  + c2%natom  
    !integer, DIMENSION(:), POINTER :: zahl
    !real (KIND=wp), DIMENSION(:), POINTER :: chrg, x, y, z
    !real (KIND=wp), DIMENSION(:), POINTER :: array
    print *, 'ADDED ', cube_add%data_title
  END FUNCTION cube_add


  FUNCTION cube_sub (c1, c2)
    TYPE(cube) :: cube_sub
    TYPE(cube), INTENT(IN) :: c1, c2
    allocate(cube_sub%array(size(c1%array)))
    cube_sub%data_title = c1%data_title // ' - ' // c2%data_title 
    cube_sub%data_desc = c1%data_desc // ' - ' // c2%data_desc 
    cube_sub%array = c1%array - c2%array
    cube_sub%natom = c1%natom  - c2%natom  
    cube_sub%nx = c1%nx
    cube_sub%ny = c1%ny
    cube_sub%nz = c1%nz
    print *, 'Subbed'
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
!  SUBROUTINE cube_cdz (mycube)
!    TYPE (cube), INTENT(IN) :: mycube
!    real, DIMENSION(size(mycube%array)) :: cube_cdz
!  END SUBROUTINE cube_cdz


  FUNCTION cube_cdz (drho)
    TYPE(cube), INTENT(IN) :: drho
    real, DIMENSION(size(drho%array)) :: cube_cdz
    real, DIMENSION(drho%nx , drho%ny , drho%nz) :: xyz
    integer :: i, x, y, z
    i = 1
    print *, 'IN cube_cdz: (', drho%nx, ', ', drho%ny, ', ', drho%nz, ')'
    !print *, 'HERE'
    do x = 1, drho%nx, 1
        do y = 1, drho%ny, 1
            do z = 1, drho%nz, 1
                xyz(x, y, z) = drho%array(i)
                i = i + 1
                if (MOD(i, 100000) .EQ. 0) then
                    print *, i, ': (', x, ', ', y, ', ', z, ') = ', drho%array(i)
                end if
            end do
        end do
    end do
    cube_cdz = drho%array
  END FUNCTION cube_cdz

END MODULE cubes
