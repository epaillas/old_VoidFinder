#!/bin/bash -f

# this script was created by vozinit, and modified to 
# correct the path and do all the work
# Rick Wagner - 15JAN08

VOBOZBIN=../bin

# build the voronoi diagrams
poe $VOBOZBIN/voz1b1_mpi single_halo.raw 0.3 1 SH_mpi 2  -nodes 1 -tasks_per_node 2
$VOBOZBIN/voztie 2 SH_mpi

# join the zones
$VOBOZBIN/jovoz adjSH_mpi.dat volSH_mpi.dat zonesSH_mpi.dat zonesSH_mpi.txt 0.01

# unbind particles
$VOBOZBIN/boz 1.0 4 0.3 1.0 single_halo.raw single_halo_vel.raw zonesSH_mpi.dat halosSH_mpi.dat halosSH_mpi.txt 1.2

cat halosSH_mpi.txt