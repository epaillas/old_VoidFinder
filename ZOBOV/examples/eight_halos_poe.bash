#!/bin/bash -f

# this script was created by vozinit, and modified to 
# correct the path and do all the work
# Rick Wagner - 15JAN08

VOBOZBIN=../bin

# build the voronoi diagrams
poe $VOBOZBIN/voz1b1_mpi eight_halos.raw 0.3 1 8H_poe 2  -nodes 1 -tasks_per_node 8
$VOBOZBIN/voztie 2 8H_poe  -nodes 1 -tasks_per_node 1

# join the zones
$VOBOZBIN/jovoz adj8H_poe.dat vol8H_poe.dat zones8H_poe.dat zones8H_poe.txt 0.01   -nodes 1 -tasks_per_node 1

# unbind particles
$VOBOZBIN/boz 1.0 0.3 1.0 eight_halos.raw eight_halos_vel.raw zones8H_poe.dat halos8H_poe.dat halos8H_poe.txt 1.2   -nodes 1 -tasks_per_node 1

cat halos8H_poe.txt