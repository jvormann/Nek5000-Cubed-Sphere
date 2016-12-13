# Nek5000-Cubed-Sphere
This Fortran program creates .rea files for the Nek5000 computational fluid dynamics code (https://nek5000.mcs.anl.gov/) that represent a cubed sphere grid (a cartesian grid is projected on a spherical shell).
To use:
-modify and run the makefile in ./src, usually you only have to set a fortran90 compiler
-run cubed.out in the main directory - for high resolutions, the files can become pretty large!
-run joinrea.sh to stitch the files together for different use scenarios (only velocity, velocity and heat, magnetohydrodynamics)
-copy the appropriate .rea file to your project folder
-remember to change the SIZE file and parameters in .rea before running
The program asks for the angular and radial resolution: the radial resolution is the number of cartesian grid points on one side of the cube, the radial resolution is the number of elements from the inner to the outer boundary. They can either be distributed equidistant or with the Gauss Lobatto distribution (denser at the boundary). Distance between the boundaries is fixed to 1.0. The total number of elements is then 6*nr*nphi*nphi, while the polynomial order and total resolution is controlled in the SIZE file of the simulation.
To create a spheroid or ellipsoid, the grid can be deformed in the usrdat2 routine in .usr.


Some untested features of the code are commented, use at your own risk!
