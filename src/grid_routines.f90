!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This module contains subroutines to create a cubed sphere grid with an
!inner shell for use with NEK5000
!Author: Jan Vormann
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE grid_routines
USE Jacobi_Subroutines
IMPLICIT NONE
CONTAINS

SUBROUTINE get_cubed_sphere (cub_sphere_coord,&
 cub_sphere_bc_T_type, cub_sphere_bc_T_value, &
&cub_sphere_bc_V_type, cub_sphere_bc_V_value, &
&cub_sphere_bc_B_type, cub_sphere_bc_B_value, &
&cub_sphere_curv,nphif, nr, r_method, phi_method,r_i,r_o, &
&mantle,dmantle,nrm,ellip, bc_t_ic, bc_t_oc,bc_t_val1, bc_t_val2)
!Returns the coordinates of the corner points of spectral elements, boundary conditions and curved sides
!in a cubed sphere grid
INTEGER :: nphif,nphi, nr, nrm
INTEGER :: r_method, phi_method, mantle ! r| 1: equidistant | 2: Gauss-Lobatto Distribution phi| 1: constant angle 2: constant distance on base cube
REAL (KIND=8), DIMENSION (1:(nr+nrm)*(nphif/4)*(nphif/4)*6,1:8,1:3) :: cub_sphere_coord !(element number, corner number, coordinate(1:x, 2:y, 3:z)
CHARACTER (len=3), DIMENSION (1:(nr+nrm)*(nphif/4)*(nphif/4)*6, 1:6) ::cub_sphere_bc_T_type
REAL (KIND=8), DIMENSION (1:(nr+nrm)*(nphif/4)*(nphif/4)*6, 1:6, 1:5) ::cub_sphere_bc_T_value
CHARACTER (len=3), DIMENSION (1:(nr+nrm)*(nphif/4)*(nphif/4)*6, 1:6) ::cub_sphere_bc_V_type
REAL (KIND=8), DIMENSION (1:(nr+nrm)*(nphif/4)*(nphif/4)*6, 1:6, 1:5) ::cub_sphere_bc_V_value
CHARACTER (len=3), DIMENSION (1:(nr+nrm)*(nphif/4)*(nphif/4)*6, 1:6) ::cub_sphere_bc_B_type
REAL (KIND=8), DIMENSION (1:(nr+nrm)*(nphif/4)*(nphif/4)*6, 1:6, 1:5) ::cub_sphere_bc_B_value
REAL (KIND=8), DIMENSION (1:(nr+nrm)*(nphif/4)*(nphif/4)*6,1:6) :: cub_sphere_curv
REAL (KIND=8), DIMENSION (6,1:(nphif/4)+1,1:(nphif/4)+1,1:3) :: cube_grid !coordinates on one side of the unit cube, x,y,z=0 in the center
REAL (KIND=8), DIMENSION (1:nr+nrm+1) :: rad_grid !complete radial coordinates for core and mantle
REAL (KIND=8), DIMENSION (1:nr+1) :: rad_grid_c !radial coordinates between r_i and r_o for the outer core
REAL (KIND=8), ALLOCATABLE, DIMENSION (:) :: rad_grid_m !radial coordinates for the mantle
REAL (KIND=8) :: r_i, r_o, l, pi !Inner and outer radius, length of side of the inscribed cube, pi
REAL (KIND=8) :: dphi, dr, veclen, dmantle,ellip
INTEGER :: in, im, ie, iv, ir, ip,ix,&
&          elcount,elcount_all,elcount_mantle,elcount_core
CHARACTER (len=3) :: bc_v_ic, bc_v_oc, bc_t_ic, bc_t_oc, bc_b_ic, bc_b_oc
CHARACTER (len=3) :: bc_v_m, bc_t_m, bc_b_m
REAL (KIND=8) :: bc_t_val1,bc_t_val2,bc_t_val3
!Set boundary condition strings
bc_v_ic='v  '
bc_v_oc='v  '
bc_v_m ='E  '
!bc_t_ic='t  '
!bc_t_oc='t  '
bc_t_m ='t  '
bc_b_ic='v  '
bc_b_oc='v  '
bc_b_m ='v  '
!bc_t_val1=0.0_8
!bc_t_val2=bc_t_val1
bc_t_val3=bc_t_val1
IF (mantle .EQ. 1) bc_b_oc='E  '
pi=4.0_8*ATAN(1.0_8)
nphi=nphif/4
dr=(r_o-r_i)/REAL(nr)
!Calculate point position on the sides of the unit cube
cube_grid=get_cube_grid(nphi, phi_method)
!Calculate the unit vector for each point
DO ie=1,6
  DO in=1,nphi+1
   DO im=1,nphi+1
    veclen=SQRT(cube_grid(ie,in,im,1)**2+cube_grid(ie,in,im,2)**2+cube_grid(ie,in,im,3)**2)
    cube_grid(ie,in,im,1:3)=cube_grid(ie,in,im,1:3)/veclen
   ENDDO
  ENDDO
ENDDO
!Calculate radial Distribution
!ALLOCATE(rad_grid(1:nr+1))
rad_grid_c=get_radial_grid(r_i, r_o, nr, r_method)
IF (mantle .EQ. 1) THEN
  ALLOCATE(rad_grid_m(1:nrm+1))
  rad_grid_m=get_radial_grid(r_o, r_o+dmantle, nrm, 3)
  rad_grid(1:nr+1)=rad_grid_c(1:nr+1)
  rad_grid(nr+2:nr+nrm+1)=rad_grid_m(2:nrm+1)
ELSE
  rad_grid(1:nr+1)=rad_grid_c(1:nr+1)
ENDIF
!Assemble the element coordinates face by face
elcount_all=6*nphi*nphi*(nr+nrm)
elcount_core=0
elcount_mantle=6*nphi*nphi*nr
IF (mantle .EQ. 1) then
        PRINT*,'Mantle count starting at: ',elcount_mantle,' of ',elcount_all
ENDIF
DO ie=1,6
 DO in=1,nphi
  DO im=1,nphi
   DO ir=1,nr+nrm
        IF (ir .LE. nr) THEN
                elcount_core=elcount_core+1
                elcount=elcount_core
        ELSE
                elcount_mantle=elcount_mantle+1
                elcount=elcount_mantle
        ENDIF
        cub_sphere_bc_T_type(elcount,1:6)='E  '
        cub_sphere_bc_T_value(elcount,1:6,1:5)=0.0_8
        cub_sphere_bc_V_type(elcount,1:6)='E  '
        cub_sphere_bc_V_value(elcount,1:6,1:5)=0.0_8
        cub_sphere_bc_B_type(elcount,1:6)='E  '
        cub_sphere_bc_B_value(elcount,1:6,1:5)=0.0_8
        !The orientation of the element depends on the face it is projected from,
        !so some IFs are needed
        IF (ie == 1) THEN
                DO ix=1,3
                        cub_sphere_coord(elcount,1,ix)=cube_grid(ie,in,im,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,2,ix)=cube_grid(ie,in+1,im,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,3,ix)=cube_grid(ie,in+1,im,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,4,ix)=cube_grid(ie,in,im,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,5,ix)=cube_grid(ie,in,im+1,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,6,ix)=cube_grid(ie,in+1,im+1,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,7,ix)=cube_grid(ie,in+1,im+1,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,8,ix)=cube_grid(ie,in,im+1,ix)*rad_grid(ir)
                ENDDO
                !Set curved faces
                cub_sphere_curv(elcount,3)=-rad_grid(ir)
                cub_sphere_curv(elcount,1)=rad_grid(ir+1)
                !Set BCs
                IF (ir == 1) THEN
                        cub_sphere_bc_T_type(elcount,3)=bc_t_ic
                        cub_sphere_bc_T_value(elcount,3,1)=bc_t_val1
                        cub_sphere_bc_V_type(elcount,3)=bc_v_ic
                        cub_sphere_bc_B_type(elcount,3)=bc_b_ic
                ELSEIF (ir == nr) THEN
                        cub_sphere_bc_T_type(elcount,1)=bc_t_oc
                        cub_sphere_bc_T_value(elcount,1,1)=bc_t_val2
                        cub_sphere_bc_V_type(elcount,1)=bc_v_oc
                        cub_sphere_bc_B_type(elcount,1)=bc_b_oc
                ELSEIF ((mantle .EQ. 1) .AND. (ir .EQ. (nr+nrm))) THEN
                        cub_sphere_bc_T_type(elcount,1)=bc_t_m
                        cub_sphere_bc_T_value(elcount,1,1)=bc_t_val3
                        cub_sphere_bc_V_type(elcount,1)=bc_v_m
                        cub_sphere_bc_B_type(elcount,1)=bc_b_m
                ENDIF
        ELSEIF (ie == 2) THEN
                DO ix=1,3
                        cub_sphere_coord(elcount,1,ix)=cube_grid(ie,in,im,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,2,ix)=cube_grid(ie,in,im,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,3,ix)=cube_grid(ie,in+1,im,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,4,ix)=cube_grid(ie,in+1,im,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,5,ix)=cube_grid(ie,in,im+1,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,6,ix)=cube_grid(ie,in,im+1,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,7,ix)=cube_grid(ie,in+1,im+1,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,8,ix)=cube_grid(ie,in+1,im+1,ix)*rad_grid(ir)
                ENDDO
                !Set curved faces
                cub_sphere_curv(elcount,4)=-rad_grid(ir)
                cub_sphere_curv(elcount,2)=rad_grid(ir+1)
                !Set BCs
                IF (ir == 1) THEN
                        cub_sphere_bc_T_type(elcount,4)=bc_t_ic
                        cub_sphere_bc_T_value(elcount,4,1)=bc_t_val1
                        cub_sphere_bc_V_type(elcount,4)=bc_v_ic
                        cub_sphere_bc_B_type(elcount,4)=bc_b_ic
                ELSEIF (ir==nr) THEN
                        cub_sphere_bc_T_type(elcount,2)=bc_t_oc
                        cub_sphere_bc_T_value(elcount,2,1)=bc_t_val2
                        cub_sphere_bc_V_type(elcount,2)=bc_v_oc
                        cub_sphere_bc_B_type(elcount,2)=bc_b_oc
                ELSEIF ((mantle .EQ. 1) .AND. (ir .EQ. (nr+nrm))) THEN
                        cub_sphere_bc_T_type(elcount,2)=bc_t_m
                        cub_sphere_bc_T_value(elcount,2,1)=bc_t_val3
                        cub_sphere_bc_V_type(elcount,2)=bc_v_m
                        cub_sphere_bc_B_type(elcount,2)=bc_b_m
                ENDIF
        ELSEIF (ie == 3) THEN
                DO ix=1,3
                        cub_sphere_coord(elcount,1,ix)=cube_grid(ie,in,im,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,2,ix)=cube_grid(ie,in+1,im,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,3,ix)=cube_grid(ie,in+1,im,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,4,ix)=cube_grid(ie,in,im,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,5,ix)=cube_grid(ie,in,im+1,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,6,ix)=cube_grid(ie,in+1,im+1,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,7,ix)=cube_grid(ie,in+1,im+1,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,8,ix)=cube_grid(ie,in,im+1,ix)*rad_grid(ir+1)
                ENDDO
                !Set curved faces
                cub_sphere_curv(elcount,1)=-rad_grid(ir)
                cub_sphere_curv(elcount,3)=rad_grid(ir+1)
                !Set BCs
                IF (ir == 1) THEN
                        cub_sphere_bc_T_type(elcount,1)=bc_t_ic
                        cub_sphere_bc_T_value(elcount,1,1)=bc_t_val1
                        cub_sphere_bc_V_type(elcount,1)=bc_v_ic
                        cub_sphere_bc_B_type(elcount,1)=bc_b_ic
                ELSEIF (ir==nr) THEN
                        cub_sphere_bc_T_type(elcount,3)=bc_t_oc
                        cub_sphere_bc_T_value(elcount,3,1)=bc_t_val2
                        cub_sphere_bc_V_type(elcount,3)=bc_v_oc
                        cub_sphere_bc_B_type(elcount,3)=bc_b_oc
                ELSEIF ((mantle .EQ. 1) .AND. (ir .EQ. (nr+nrm))) THEN
                        cub_sphere_bc_T_type(elcount,3)=bc_t_m
                        cub_sphere_bc_T_value(elcount,3,1)=bc_t_val3
                        cub_sphere_bc_V_type(elcount,3)=bc_v_m
                        cub_sphere_bc_B_type(elcount,3)=bc_b_m
                ENDIF
        ELSEIF (ie == 4) THEN
                DO ix=1,3
                        cub_sphere_coord(elcount,1,ix)=cube_grid(ie,in,im,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,2,ix)=cube_grid(ie,in,im,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,3,ix)=cube_grid(ie,in+1,im,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,4,ix)=cube_grid(ie,in+1,im,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,5,ix)=cube_grid(ie,in,im+1,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,6,ix)=cube_grid(ie,in,im+1,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,7,ix)=cube_grid(ie,in+1,im+1,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,8,ix)=cube_grid(ie,in+1,im+1,ix)*rad_grid(ir+1)
                ENDDO
                !Set curved faces
                cub_sphere_curv(elcount,2)=-rad_grid(ir)
                cub_sphere_curv(elcount,4)=rad_grid(ir+1)
                !Set BCs
                IF (ir == 1) THEN
                        cub_sphere_bc_T_type(elcount,2)=bc_t_ic
                        cub_sphere_bc_T_value(elcount,2,1)=bc_t_val1
                        cub_sphere_bc_V_type(elcount,2)=bc_v_ic
                        cub_sphere_bc_B_type(elcount,2)=bc_b_ic
                ELSEIF (ir==nr) THEN
                        cub_sphere_bc_T_type(elcount,4)=bc_t_oc
                        cub_sphere_bc_T_value(elcount,4,1)=bc_t_val2
                        cub_sphere_bc_V_type(elcount,4)=bc_v_oc
                        cub_sphere_bc_B_type(elcount,4)=bc_b_oc
                ELSEIF ((mantle .EQ. 1) .AND. (ir .EQ. (nr+nrm))) THEN
                        cub_sphere_bc_T_type(elcount,4)=bc_t_m
                        cub_sphere_bc_T_value(elcount,4,1)=bc_t_val3
                        cub_sphere_bc_V_type(elcount,4)=bc_v_m
                        cub_sphere_bc_B_type(elcount,4)=bc_b_m
                ENDIF
        ELSEIF (ie == 5) THEN
                DO ix=1,3
                   cub_sphere_coord(elcount,1,ix)=cube_grid(ie,in,im,ix)*rad_grid(ir+1)
                   cub_sphere_coord(elcount,2,ix)=cube_grid(ie,in+1,im,ix)*rad_grid(ir+1)
                   cub_sphere_coord(elcount,3,ix)=cube_grid(ie,in+1,im+1,ix)*rad_grid(ir+1)
                   cub_sphere_coord(elcount,4,ix)=cube_grid(ie,in,im+1,ix)*rad_grid(ir+1)
                   cub_sphere_coord(elcount,5,ix)=cube_grid(ie,in,im,ix)*rad_grid(ir)
                   cub_sphere_coord(elcount,6,ix)=cube_grid(ie,in+1,im,ix)*rad_grid(ir)
                   cub_sphere_coord(elcount,7,ix)=cube_grid(ie,in+1,im+1,ix)*rad_grid(ir)
                   cub_sphere_coord(elcount,8,ix)=cube_grid(ie,in,im+1,ix)*rad_grid(ir)
                ENDDO
                !Set curved faces
                cub_sphere_curv(elcount,6)=-rad_grid(ir)
                cub_sphere_curv(elcount,5)=rad_grid(ir+1)
                !Set BCs
                IF (ir == 1) THEN
                        cub_sphere_bc_T_type(elcount,6)=bc_t_ic
                        cub_sphere_bc_T_value(elcount,6,1)=bc_t_val1
                        cub_sphere_bc_V_type(elcount,6)=bc_v_ic
                        cub_sphere_bc_B_type(elcount,6)=bc_b_ic
                ELSEIF (ir==nr) THEN
                        cub_sphere_bc_T_type(elcount,5)=bc_t_oc
                        cub_sphere_bc_T_value(elcount,5,1)=bc_t_val2
                        cub_sphere_bc_V_type(elcount,5)=bc_v_oc
                        cub_sphere_bc_B_type(elcount,5)=bc_b_oc
                ELSEIF ((mantle .EQ. 1) .AND. (ir .EQ. (nr+nrm))) THEN
                        cub_sphere_bc_T_type(elcount,5)=bc_t_m
                        cub_sphere_bc_T_value(elcount,5,1)=bc_t_val3
                        cub_sphere_bc_V_type(elcount,5)=bc_v_m
                        cub_sphere_bc_B_type(elcount,5)=bc_b_m
                ENDIF
        ELSEIF (ie == 6) THEN
                DO ix=1,3
                        cub_sphere_coord(elcount,1,ix)=cube_grid(ie,in,im,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,2,ix)=cube_grid(ie,in+1,im,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,3,ix)=cube_grid(ie,in+1,im+1,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,4,ix)=cube_grid(ie,in,im+1,ix)*rad_grid(ir)
                        cub_sphere_coord(elcount,5,ix)=cube_grid(ie,in,im,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,6,ix)=cube_grid(ie,in+1,im,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,7,ix)=cube_grid(ie,in+1,im+1,ix)*rad_grid(ir+1)
                        cub_sphere_coord(elcount,8,ix)=cube_grid(ie,in,im+1,ix)*rad_grid(ir+1)
                ENDDO
                !Set curved faces
                cub_sphere_curv(elcount,5)=-rad_grid(ir)
                cub_sphere_curv(elcount,6)=rad_grid(ir+1)
                !Set BCs
                IF (ir == 1) THEN
                        cub_sphere_bc_T_type(elcount,5)=bc_t_ic
                        cub_sphere_bc_T_value(elcount,5,1)=bc_t_val1
                        cub_sphere_bc_V_type(elcount,5)=bc_v_ic
                        cub_sphere_bc_B_type(elcount,5)=bc_b_ic
                ELSEIF (ir==nr) THEN
                        cub_sphere_bc_T_type(elcount,6)=bc_t_oc
                        cub_sphere_bc_T_value(elcount,6,1)=bc_t_val2
                        cub_sphere_bc_V_type(elcount,6)=bc_v_oc
                        cub_sphere_bc_B_type(elcount,6)=bc_b_oc
                ELSEIF ((mantle .EQ. 1) .AND. (ir .EQ. (nr+nrm))) THEN
                        cub_sphere_bc_T_type(elcount,6)=bc_t_m
                        cub_sphere_bc_T_value(elcount,6,1)=bc_t_val3
                        cub_sphere_bc_V_type(elcount,6)=bc_v_m
                        cub_sphere_bc_B_type(elcount,6)=bc_b_m
                ENDIF
          ENDIF
!         Set ellipticity
!          cub_sphere_coord(:,:,3)=ellip*cub_sphere_coord(:,:,3)
!         Fix all element sides in the mantle for additional stability
!          IF ((ir .gt. nr) .AND. (ir .LE. (nr+nrm))) THEN
!                cub_sphere_bc_V_type(elcount,1:6)=bc_v_m
!          ENDIF
    ENDDO
  ENDDO
 ENDDO
ENDDO
PRINT*,'Element # core: ',elcount_core
IF (mantle .EQ. 1) PRINT*,'Element # mantle: ',elcount_mantle
RETURN
END SUBROUTINE get_cubed_sphere

FUNCTION get_cube_grid(nphi, phi_method)
!Returns the distribution of points at the faces of a unit cube
!phi_method=1 -> Points have the same angular distance when projected on a circle
!phi_method=2 -> Points have the same linear distance parallel to the edges
REAL (KIND=8) :: l, pi, dphi
REAL (KIND=8), DIMENSION(6,1:nphi+1,1:nphi+1,1:3) :: get_cube_grid
REAL (KIND=8), DIMENSION(1:nphi+1) :: grid_1d
INTEGER :: phi_method
INTEGER :: nphi, in, im
pi=4.0_8*ATAN(1.0_8)
l=1.0_8
!Start by creating a 1D Distribution
IF (phi_method == 1) THEN
dphi=pi/(2.0_8*REAL(nphi))
  IF (MOD(nphi, 2) .NE. 0.) THEN
    DO in=((nphi+1)/2)+1, nphi+1
      IF (in==((nphi+1)/2)+1) THEN
        grid_1d(in)=l/2.0_8*TAN(0.5_8*dphi)
      ELSEIF (in==((nphi+1)/2)+2) THEN
        grid_1d(in)=l/2.0_8*TAN(1.5_8*dphi)
      ELSE
        grid_1d(in)=l/2.0_8*TAN(0.5_8*dphi+in-(REAL((nphi+1)/2)-1.0_8)*dphi)
      ENDIF
       grid_1d(nphi+1-(in-1))=-grid_1d(in)
    ENDDO
  ELSE
    DO in=(nphi+1)/2,nphi+1
      IF (in==(nphi+1)/2) THEN
        grid_1d(in)=0.0_8
      ELSE
        grid_1d(in)=l/2.0_8*TAN((in-(nphi+1)/2-1)*dphi)
        grid_1d((nphi+1)-(in-1))=-grid_1d(in)
      ENDIF
    ENDDO
  ENDIF
 ELSEIF (phi_method == 2) THEN
 dphi=l/nphi
 DO in=1,nphi+1
  grid_1d(in)=-l/2.0_8+REAL(in-1)*dphi
  !PRINT*,grid_1d(in)
 ENDDO
ELSE
 PRINT*,'Method for angular grid invalid.'
 get_cube_grid=0.0_8
 STOP
ENDIF
!Set constant coordinates first
!Face 1:Front, Face 3: Back
get_cube_grid(1,:,:,2)=-0.5_8
get_cube_grid(3,:,:,2)=+0.5_8
!Face 4:Left, Face 2: Right
get_cube_grid(4,:,:,1)=-0.5_8
get_cube_grid(2,:,:,1)=+0.5_8
!Face 6:Top, Face 5: Bottom
get_cube_grid(6,:,:,3)=+0.5_8
get_cube_grid(5,:,:,3)=-0.5_8
!Set variable coordinates
DO in=1,nphi+1
  DO im=1,nphi+1
   get_cube_grid(1,in,im,1)=grid_1d(in)
   get_cube_grid(1,in,im,3)=grid_1d(im)
   get_cube_grid(6,in,im,1)=grid_1d(in)
   get_cube_grid(6,in,im,2)=grid_1d(im)
   get_cube_grid(3,in,im,1)=grid_1d(in)
   get_cube_grid(3,in,im,3)=grid_1d(im)
   get_cube_grid(4,in,im,2)=grid_1d(in)
   get_cube_grid(4,in,im,3)=grid_1d(im)
   get_cube_grid(2,in,im,2)=grid_1d(in)
   get_cube_grid(2,in,im,3)=grid_1d(im)
   get_cube_grid(5,in,im,1)=grid_1d(in)
   get_cube_grid(5,in,im,2)=grid_1d(im)
  ENDDO
ENDDO
END FUNCTION get_cube_grid
!
!
FUNCTION get_radial_grid(ri, ro, nr, r_method)
REAL (KIND=8), DIMENSION (1:nr+1) :: get_radial_grid,gll_qp
REAL (KIND=8), DIMENSION (1:nr-1) :: gll_qp_base
REAL (KIND=8) :: ri, ro, dr, asum, a
INTEGER :: nr, r_method, ir
IF (r_method==1) THEN
  dr=(ro-ri)/REAL(nr)
  DO ir=1,nr+1
    get_radial_grid(ir)=ri+REAL(ir-1)*dr
  ENDDO
ELSEIF (r_method==2) THEN
    gll_qp=REAL(gll_quad_points(nr+1))
    get_radial_grid=map_ltog(nr+1, gll_qp, ri, ro)
ELSEIF (r_method==3) THEN !Only used for the mantle grid, which becomes coarser with increasing radius.
    a=1.50_8
    asum=0.0_8
    DO ir=1,nr
      asum=asum+a**(ir-1)
    ENDDO
    get_radial_grid(1)=ri
    get_radial_grid(2)=(ro+ri*(asum-1))/asum
    dr=get_radial_grid(2)-get_radial_grid(1)
    DO ir=3,nr
     get_radial_grid(ir)=get_radial_grid(ir-1)+a**(ir-2)*dr
    ENDDO
    get_radial_grid(nr+1)=ro
    DO ir=1,nr+1
     PRINT*, ir,get_radial_grid(ir)
    ENDDO
ELSE
 PRINT*,'Method for radial grid invalid.'
 get_radial_grid=0.0_8
 STOP
ENDIF
RETURN
END FUNCTION get_radial_grid
!
FUNCTION map_ltog(q,gll_qp,xl,xr)
!This function maps the local coordinates gll_qp, given on the standard element [-1:1] to any
!global element [xl:xr]. gll_qp should be the Gauss-Lobatto-Legendre-Quadrature Points
INTEGER q, iq
REAL (KIND=8), DIMENSION (0:q-1) :: gll_qp, map_ltog
REAL (KIND=8) :: xl, xr
DO iq=0,q-1
  map_ltog(iq)=(1.0_8-gll_qp(iq))/2.0_8*xl+(1.0_8+gll_qp(iq))/2.0_8*xr
ENDDO
RETURN
END FUNCTION
!
SUBROUTINE write_elements_to_rea(grid,nel,nelcore,nelmantle,element_file)
!Writes the element definition part of the .rea file to element_file,
!using the elements defined in grid (only for 3D)
INTEGER :: nel,nelcore,nelmantle,iel
INTEGER :: ascii_count, n_count
REAL (KIND=8), DIMENSION (1:nel,1:8,1:3) :: grid
CHARACTER (len=100) :: element_file
OPEN(10,file=element_file)
WRITE(10,'(A)') ' **MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8'
WRITE(10,'(3I8,1X,A)') nel,3,nelcore,'NEL,NDIM,NELV'
ascii_count=97
n_count=1
DO iel=1,nel
 WRITE(10,'(18X,I6,4X,I3,A1,11x,i5)') iel,n_count,ACHAR(ascii_count),0
 WRITE(10,*) grid(iel,1,1),grid(iel,2,1),grid(iel,3,1),grid(iel,4,1)
 WRITE(10,*) grid(iel,1,2),grid(iel,2,2),grid(iel,3,2),grid(iel,4,2)
 WRITE(10,*) grid(iel,1,3),grid(iel,2,3),grid(iel,3,3),grid(iel,4,3)
 WRITE(10,*) grid(iel,5,1),grid(iel,6,1),grid(iel,7,1),grid(iel,8,1)
 WRITE(10,*) grid(iel,5,2),grid(iel,6,2),grid(iel,7,2),grid(iel,8,2)
 WRITE(10,*) grid(iel,5,3),grid(iel,6,3),grid(iel,7,3),grid(iel,8,3)
 ascii_count=ascii_count+1
 IF (ascii_count == 123) THEN
  ascii_count=97
  n_count=n_count+1
 ENDIF
ENDDO
CLOSE(10)
END SUBROUTINE
!
SUBROUTINE write_curvature_to_rea(curves,nel,curvature_file)
!Writes the curvature definition part of the .rea file to curvature_file,
!using the elements defined in grid (only for 3D)
INTEGER :: nel, iel, iface
REAL (KIND=8), DIMENSION (1:nel,1:6) :: curves
CHARACTER (len=100) :: curvature_file, my_format

OPEN(10,file=curvature_file)
WRITE(10,*) '  ***** CURVED SIDE DATA *****'
WRITE(10,'(I8,A)') nel*2,' Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE'
IF (nel .LT. 1000) THEN
  my_format='(I3,I3,1p5G14.6, 1X,A1)'
ELSEIF (nel .LT. 1000000) THEN
  my_format='(I2,I6,1p5G14.6, 1X,A1)'
ELSE
  my_format='(I2,I10,1p5G14.6, 1X,A1)'
ENDIF
DO iel=1,nel
 DO iface=1,6
  IF (curves(iel,iface) .NE. 0.0_8) THEN
    WRITE(10,my_format) iface,iel,0.,0.,0.,ABS(curves(iel,iface)),1.,'s'
  ENDIF
 ENDDO
ENDDO
CLOSE(10)
END SUBROUTINE
!
SUBROUTINE write_bcs_to_rea(bc_type,bc_value,nel,nelcore,nelmantle,bc_file,bc_var)
INTEGER :: nel,nelcore,nelmantle, iel, iface
INTEGER :: nelidx
character(*), DIMENSION (1:nel,1:6) :: bc_type
REAL (KIND=8), DIMENSION (1:nel, 1:6, 1:5) :: bc_value
CHARACTER (len=100) :: bc_file, my_format
CHARACTER(*) :: bc_var

OPEN(10,file=bc_file)
IF (bc_var == 'T') THEN
  WRITE(10,'(A)') '  ***** THERMAL BOUNDARY CONDITIONS *****'
ELSEIF (bc_var == 'V') THEN
  WRITE(10,'(A)') '  ***** BOUNDARY CONDITIONS *****'
  WRITE(10,'(A)') '  ***** FLUID   BOUNDARY CONDITIONS *****'
ELSEIF (bc_var == 'B') THEN
  WRITE(10,'(A)') '  ***** PASSIVE SCALAR 1  BOUNDARY CONDITIONS *****'
ELSE
  PRINT*,'BC Type ',bc_var,' not valid in write_bcs_to_rea'
  STOP
ENDIF
IF (nel .LT. 1000) THEN
  my_format='(1X,A3,2I3,5G14.6)'
ELSEIF (nel .LT. 100000) THEN
  my_format='(1X,A3,I5,I1,5G14.6)'
ELSE
  my_format='(1X,A3,I10,I1,5G14.6)'
ENDIF
nelidx=nel
IF (bc_var == 'V') nelidx=nelcore
DO iel=1,nelidx
        DO iface=1,6
            WRITE(10,my_format) bc_type(iel,iface),iel,iface,bc_value(iel,iface,1:5)
        ENDDO
ENDDO
END SUBROUTINE

!
END MODULE
