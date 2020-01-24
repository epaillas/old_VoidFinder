#!/bin/bash -f

# this script was created by vozinit, and modified to 
# correct the path and do all the work
# Rick Wagner - 15JAN08

VOBOZBIN=../bin

# build the voronoi diagrams
mpirun -np 2 $VOBOZBIN/voz1b1_mpi single_halo_boundary.raw 0.3 1 SHB_mpi 2  
$VOBOZBIN/voztie 2 SHB_mpi

# join the zones
$VOBOZBIN/jovoz adjSHB_mpi.dat volSHB_mpi.dat zonesSHB_mpi.dat zonesSHB_mpi.txt 0.01

# unbind particles
$VOBOZBIN/boz 1.0 4 0.3 1.0 single_halo_boundary.raw single_halo_boundary_vel.raw zonesSHB_mpi.dat halosSHB_mpi.dat halosSHB_mpi.txt 1.2

cat halosSHB_mpi.txt