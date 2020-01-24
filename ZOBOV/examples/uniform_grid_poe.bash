#!/bin/bash -f

# this script was created by vozinit, and modified to 
# correct the path and do all the work
# Rick Wagner - 15JAN08

VOBOZBIN=../bin

# build the voronoi diagrams
poe $VOBOZBIN/voz1b1_mpi uniform_grid.raw 0.3 1 UG_mpi 2 -nodes 1 -tasks_per_node 2
$VOBOZBIN/voztie 2 UG_mpi

# join the zones
$VOBOZBIN/jovoz adjUG_mpi.dat volUG_mpi.dat zonesUG_mpi.dat zonesUG_mpi.txt 0.01

# unbind particles
$VOBOZBIN/boz 1.0 0.3 1.0 uniform_grid.raw uniform_grid_vel.raw zonesUG_mpi.dat halosUG_mpi.dat halosUG_mpi.txt 1.2

cat halosUG_mpi.txt