#!/bin/bash -f

# this script was created by vozinit, and modified to 
# correct the path and do all the work
# Rick Wagner - 15JAN08

VOBOZBIN=../bin

# build the voronoi diagrams
mpirun -np 2 $VOBOZBIN/voz1b1_mpi eight_halos.raw 0.3 1 8H_mpi 2
$VOBOZBIN/voztie 2 8H_mpi

# join the zones
$VOBOZBIN/jovoz adj8H_mpi.dat vol8H_mpi.dat zones8H_mpi.dat zones8H_mpi.txt 0.01

# unbind particles
$VOBOZBIN/boz 1.0 0.3 1.0 eight_halos.raw eight_halos_vel.raw zones8H_mpi.dat halos8H_mpi.dat halos8H_mpi.txt 1.2

cat halos8H_mpi.txt