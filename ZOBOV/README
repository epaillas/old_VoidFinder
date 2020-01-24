= About =

This package includes VOBOZ (VOronoi BOund Zones), version 1.3.2, a
halo-finding algorithm by Mark Neyrinck (with "advice" from
grad-school advisors Andrew Hamilton and Nick Gnedin).  It also
includes the void-finding algorithm ZOBOV.  Many thanks to Rick Wagner
(http://lca.ucsd.edu/projects/rpwagner/; rick@ucsd.edu) for doing many
modifications to bring VOBOZ to v.1.2 and 1.3.  Some of this README
file is plagiarized from him, too.  Documentation on VOBOZ is at

  http://www.ifa.hawaii.edu/~neyrinck/voboz.  

Rick also put up a nice wiki containing some information about the
examples he created, and how to run VOBOZ:

  http://lca.ucsd.edu/projects/rpwagner/wiki/VOBOZ

This version contains both the VOBOZ/ZOBOV source code, and the Qhull code.

For help in understanding VOBOZ/ZOBOV, also see the papers describing them, at
http://arxiv.org/abs/astro-ph/0402346 (VOBOZ; Neyrinck, Gnedin & Hamilton), and
http://arxiv.org/abs/0712.3049 (ZOBOV; Neyrinck).
If these resources are insufficient, feel free to email
Mark.Neyrinck@colorado.edu with questions.

See versions.txt for a version history.

= Contents =

qhull    - Qhull source code
	   Qhull has not been modified, and is included only for convenience.
voboz    - VOBOZ source code
bin      - Contains some python scripts; also, 
           VOBOZ binaries get installed there
examples - code to generate sample datasets and scripts

= Building = 

To build VOBOZ, modify the Makefiles in src/, qhull/src/ and examples/ for 
your local machine. Currently, the files are set up to use GCC on generic
Unix/Linux type operating systems (including OS X).

To build, run `make` from this directory. This will compile qhull and the
VOBOZ programs, and place the VOBOZ programs in bin/. `make clean` removes
all binary files, etc.

= Examples =

Running `make examples` will build the sample data in examples/, and run 
VOBOZ on it. Currently, these examples are mostly for reference on how to run
the programs, and how data is handled.

= Precision =

The default precision is double (8 byte floats); to change this, edit 
qhull/src/user.h, and change REALfloat to 1.

= Mini-License =

This is free software.  It may be freely copied, modified, and
redistributed, as long as the authors are acknowledged.  There is no
warranty or other guarantee of fitness for VOBOZ or ZOBOV; they are
provided "as are."  We welcome bug reports (but do not guarantee their
fixing), which should be sent to Mark.Neyrinck@colorado.edu.

File last modified Mar 21 2009
