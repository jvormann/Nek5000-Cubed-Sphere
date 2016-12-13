PROGRAM cubed_sphere
USE grid_routines
IMPLICIT NONE
INTEGER :: nphi, nr, nel, nrm, nelcore, nelmantle
INTEGER :: r_method, phi_method, mantle
REAL (KIND=8), ALLOCATABLE, DIMENSION (:,:,:) :: grid
CHARACTER (len=3), ALLOCATABLE, DIMENSION (:,:) :: grid_bc_T_type,  grid_bc_V_type, grid_bc_B_type
REAL (KIND=8), ALLOCATABLE, DIMENSION (:,:,:) :: grid_bc_T_value, grid_bc_V_value, grid_bc_B_value
Real (KIND=8), ALLOCATABLE, DIMENSION (:,:) :: grid_curv
REAL (KIND=8) :: r_i, r_o, d, shell_ratio, dmantle
REAL (KIND=8) :: ellip
CHARACTER (len=100) :: element_file, curvature_file, bc_T_file, bc_V_file, bc_B_file
PRINT*,'Number of angular elements NPHI for one side of the cube?'
READ*,nphi
nphi=4*nphi
PRINT*,'Number of radial elements NR?'
READ*,nr
IF ( nr .LE. 1) nr=2
!Option (2) uses an equidistant grid on the cube sides. (1) creates a grid with equal
!angles when projecting. Usually, (1) is needed. Uncomment to enable option.
!PRINT*,'Angular Distribution: constant angle (1) or constant linear distance (2)'
!READ*,phi_method
phi_method=2
IF ((phi_method .NE. 1) .AND. (phi_method .NE. 2)) STOP
PRINT*,'Radial Distribution: equidistant (1) or Gauss Lobatto (2)'
READ*, r_method
IF ((r_method .NE. 1) .AND. (r_method .NE. 2)) STOP
PRINT*,'Ratio inner/outer radius?'
READ*,shell_ratio
!The thickness of the shell should usually be 1. Uncomment to enable option.
!PRINT*,'Spherical shell thickness?'
!READ*,d
d=1.0_8
!Changes to the ellipticity should be made in the .usr file (usrdat2 routine)
!Uncomment to enable option.
!PRINT*,'Ellipticity c/a?'
!READ*,ellip
ellip=1.0d0
r_o=d/(1.0_8-shell_ratio)
r_i=r_o-d
PRINT*,'Outer radius r_o=',r_o
PRINT*,'Inner radius r_i=',r_i
!This add elements to use for the thermal field only.
!Uncomment to enable option.
!PRINT*,'Add a mantle layer? (1) Yes -- (2) No'
!READ*,mantle
nrm=0
mantle=2
IF (mantle .EQ. 1) THEN
  PRINT*,'NR for the mantle?'
  READ*,nrm
  PRINT*,'Mantle thickness?'
  READ*,dmantle
ENDIF
nel=(nr+nrm)*(nphi/4)*(nphi/4)*6
nelcore=nr*(nphi/4)*(nphi/4)*6
nelmantle=nrm*(nphi/4)*(nphi/4)*6
PRINT*,'Total Elements: ', nel
IF (mantle .EQ. 1) THEN
        PRINT*,'Core Elements: ', nelcore
        PRINT*,'Mantle Elements: ', nelmantle
        PRINT*,'Remember to change *.rea accordingly!'
ENDIF
ALLOCATE(grid(1:nel,1:8,1:3))
ALLOCATE(grid_bc_T_type(nel,1:6),grid_bc_V_type(nel,1:6),grid_bc_B_type(nel,1:6))
ALLOCATE(grid_bc_T_value(nel,1:6,1:5),grid_bc_V_value(nel,1:6,1:5),grid_bc_B_value(nel,1:6,1:5))
ALLOCATE(grid_curv(1:nel,1:6))
CALL get_cubed_sphere(grid,grid_bc_T_type, grid_bc_T_value,grid_bc_V_type, &
                          &grid_bc_V_value,grid_bc_B_type, grid_bc_B_value, &
                          &grid_curv,nphi, nr, r_method, phi_method, r_i, r_o, &
                          &mantle,dmantle,nrm,ellip)
element_file='elements.rea'
curvature_file='curvature.rea'
bc_V_file='bc_V.rea'
bc_T_file='bc_T.rea'
bc_B_file='bc_B.rea'
CALL write_elements_to_rea(grid, nel, nelcore, nelmantle, element_file)
CALL write_curvature_to_rea(grid_curv, nel, curvature_file)
CALL write_bcs_to_rea(grid_bc_V_type, grid_bc_V_value, nel, nelcore, nelmantle, bc_V_file, 'V')
CALL write_bcs_to_rea(grid_bc_T_type, grid_bc_T_value, nel, nelcore, nelmantle, bc_T_file, 'T')
CALL write_bcs_to_rea(grid_bc_B_type, grid_bc_B_value, nel, nelcore, nelmantle, bc_B_file, 'B')
WRITE(*,'(A)') 'Done, call joinrea.sh to create .rea files'
END PROGRAM cubed_sphere
