#!/bin/bash -f

# this script was created by vozinit, and modified to 
# correct the path and do all the work
# Rick Wagner - 15JAN08

VOBOZBIN=../bin

# build the voronoi diagrams
poe $VOBOZBIN/voz1b1_serial eight_halos.raw 0.3 1 8H 2 0 0 0  -nodes 1 -tasks_per_node 1
poe $VOBOZBIN/voz1b1_serial eight_halos.raw 0.3 1 8H 2 0 0 1  -nodes 1 -tasks_per_node 1
poe $VOBOZBIN/voz1b1_serial eight_halos.raw 0.3 1 8H 2 0 1 0  -nodes 1 -tasks_per_node 1
poe $VOBOZBIN/voz1b1_serial eight_halos.raw 0.3 1 8H 2 0 1 1  -nodes 1 -tasks_per_node 1
poe $VOBOZBIN/voz1b1_serial eight_halos.raw 0.3 1 8H 2 1 0 0  -nodes 1 -tasks_per_node 1
poe $VOBOZBIN/voz1b1_serial eight_halos.raw 0.3 1 8H 2 1 0 1  -nodes 1 -tasks_per_node 1
poe $VOBOZBIN/voz1b1_serial eight_halos.raw 0.3 1 8H 2 1 1 0  -nodes 1 -tasks_per_node 1
poe $VOBOZBIN/voz1b1_serial eight_halos.raw 0.3 1 8H 2 1 1 1  -nodes 1 -tasks_per_node 1
poe $VOBOZBIN/voztie 2 8H

# join the zones
poe $VOBOZBIN/jovoz adj8H.dat vol8H.dat zones8H.dat zones8H.txt 0.01  -nodes 1 -tasks_per_node 1

# unbind particles
poe $VOBOZBIN/boz 1.0 0.3 1.0 eight_halos.raw eight_halos_vel.raw zones8H.dat halos8H.dat halos8H.txt 1.2  -nodes 1 -tasks_per_node 1

cat halos8H.txt
