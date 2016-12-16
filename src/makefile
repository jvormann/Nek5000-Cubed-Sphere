FCC = gfortran
#FCC = xlf
FILES = jacobi.f90 grid_routines.f90 cubed_sphere.f90
OUTFILE = cubed.out
all: $(FILES)
	$(FCC)  -o ../$(OUTFILE)  -static $(FILES)

clean:
	rm -f grid_routines.mod jacobi_subroutines.mod ../$(OUTFILE)
