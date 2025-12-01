PROGRAM cda
  USE kinds, ONLY: wp => dp
  USE cubes
  IMPLICIT NONE

  ! 1. Insantiate cube objects with data from cube files
  ! ../test/CuCO+/ab.cube, ../test/CuCO+/a.cube, and ../test/CuCO+/b.cube
  TYPE(cube) :: rho_a, rho_b, rho_ab, rho_ref, drho
  REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: cdz
  call cube_get(rho_a, '../test/CuCO+/a.cube')
  call cube_get(rho_b, '../test/CuCO+/b.cube')
  call cube_get(rho_ab, '../test/CuCO+/ab.cube')
  !call cube_print(rho_a)
  !call cube_print(rho_b)

  ! 2. Generate a cube object containing the charge redistribution:
  ! \Delta \rho = \rho_{AB} - \rho_{A} - \rho_{B}
  print *, 'Adding rho_a & rho_b . . .'
  rho_ref = rho_a + rho_b
  call cube_print(rho_ref)
  print *, 'Subtracting rho_ab - rho_ref . . .'
  drho = rho_ab - rho_ref
  call cube_print(drho)

  ! 3. Use an ad hoc defined external procedures in module cubes for extracting
  ! the charge-displacement function from the charge redistribution
  cdz = cube_cdz(drho)

END PROGRAM cda
