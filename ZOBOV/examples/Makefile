# Compiler choice:

# Gcc
CC = gcc
CFLAGS = -g -O2

QINC = -I../qhull/src

###############

EXES = uniform_grid single_halo eight_halos single_halo_boundary
SCRIPTS =  $(patsubst %,%.bash,$(EXES)) $(patsubst %,%_mpi.bash,$(EXES)) 

all: $(EXES)

data: particles halos

uniform_grid: uniform_grid.o
	$(CC) $(CFLAGS) -o $@ $@.o

single_halo: single_halo.o
	$(CC) $(CFLAGS) -o $@ $@.o

single_halo_boundary: single_halo_boundary.o
	$(CC) $(CFLAGS) -o $@ $@.o

eight_halos: eight_halos.o
	$(CC) $(CFLAGS) -o $@ $@.o

particles: $(EXES)
	for EXE in $(EXES); \
	do ./$$EXE; done	

halos: particles
	for SCRIPT in $(SCRIPTS); \
	do ./$$SCRIPT; done

.c.o:
	$(CC) -c -o $@ $(CFLAGS) $(QINC) $*.c

clean:
	rm -f $(EXES) *.o *~ *.raw *.dat *.txt part.*

