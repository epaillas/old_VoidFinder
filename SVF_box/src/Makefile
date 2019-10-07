# Compiler choice

FC = gfortran
MPIFC = mpif90

MLIBS = -lm

#################

EXES_3D = grow_spheres.exe recentring.exe overlapping.exe vg_ccf_monopole.exe
EXES_2D = grow_spheres_2D.exe recentring_2D.exe overlapping_2D.exe
EXES = $(EXES_3D) $(EXES_2D)

all: $(EXES)

grow_spheres.exe: grow_spheres.o
	$(MPIFC) grow_spheres.o -o grow_spheres.exe $(MLIBS)

recentring.exe: recentring.o
	$(MPIFC) recentring.o -o recentring.exe $(MLIBS)

overlapping.exe: overlapping.o
	$(MPIFC) overlapping.o -o overlapping.exe $(MLIBS)

vg_ccf_monopole.exe: vg_ccf_monopole.o
	$(MPIFC) vg_ccf_monopole.o -o vg_ccf_monopole.exe $(MLIBS)

grow_spheres.o: grow_spheres.f90
	$(MPIFC) -c grow_spheres.f90

recentring.o: recentring.f90
	$(MPIFC) -c recentring.f90

overlapping.o: overlapping.f90
	$(MPIFC) -c overlapping.f90

vg_ccf_monopole.o: vg_ccf_monopole.f90
	$(MPIFC) -c vg_ccf_monopole.f90


grow_spheres_2D.exe: grow_spheres_2D.o
	$(MPIFC) grow_spheres_2D.o -o grow_spheres_2D.exe $(MLIBS)

recentring_2D.exe: recentring_2D.o
	$(MPIFC) recentring_2D.o -o recentring_2D.exe $(MLIBS)

overlapping_2D.exe: overlapping_2D.o
	$(MPIFC) overlapping_2D.o -o overlapping_2D.exe $(MLIBS)

grow_spheres_2D.o: grow_spheres_2D.f90
	$(MPIFC) -c grow_spheres_2D.f90

recentring_2D.o: recentring_2D.f90
	$(MPIFC) -c recentring_2D.f90

overlapping_2D.o: overlapping_2D.f90
	$(MPIFC) -c overlapping_2D.f90

r8lib.o: r8lib.f90
	$(FC) -c r8lib.f90

install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod
