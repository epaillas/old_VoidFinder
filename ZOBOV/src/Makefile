# Compiler choice:

# Gcc
CC    = gcc
MPICC = mpicc
CFLAGS = -O2 -g -ansi

MLIBS	=   -lm

QLIB = -L../qhull/src/ -lqhull
QINC = -I../qhull/src/

###############

EXES = vozinit voz1b1 voztie jovoz boz jozov vozisol ascii_to_bin
EXES_NOMPI = vozinit voz1b1 voztie jovoz boz jozov vozisol 

all: $(EXES)

all_nompi: $(EXES_NOMPI)

jovoz: jovoz.o findrtop.o voz.h
	$(CC) $(CFLAGS) -o jovoz jovoz.o findrtop.o $(MLIBS)

jozov: jozov.o findrtop.o
	$(CC) $(CFLAGS) -o jozov jozov.o findrtop.o $(MLIBS)

boz: boz.o readfiles.o findrtop.o voz.h
	$(CC) $(CFLAGS) -o boz boz.o readfiles.o findrtop.o $(MLIBS)

voz1b1: voz1b1_serial.o voz1b1.o readfiles.o vozutil.o  voz.h
	$(CC)  -o voz1b1 $(CFLAGS) voz1b1_serial.o voz1b1.o readfiles.o vozutil.o -L. $(QLIB) $(MLIBS)

vozinit: vozinit.o readfiles.o voz.h
	$(CC)  -o vozinit $(CFLAGS) vozinit.o readfiles.o -L. $(MLIBS) 

voztie: voztie.o readfiles.o voz.h
	$(CC)  -o voztie $(CFLAGS) voztie.o readfiles.o

vozisol: vozisol.o readfiles.o vozutil.o  voz.h
	$(CC)  -o vozisol $(CFLAGS) vozisol.o readfiles.o vozutil.o -L. $(QLIB) $(MLIBS)

ascii_to_bin: ascii_to_bin.o 
	$(CC) $(CFLAGS) -o ascii_to_bin ascii_to_bin.o $(MLIBS)

install: $(EXES)
	cp $(EXES) ../bin

.c.o:
	$(CC) -c -o $@ $(CFLAGS) $(QINC) $*.c

clean:
	rm -f $(EXES) *.o *~
