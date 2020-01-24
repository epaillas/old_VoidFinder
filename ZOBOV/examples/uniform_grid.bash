#!/bin/bash -f

# this script was created by vozinit, and modified to 
# correct the path and do all the work
# Rick Wagner - 15JAN08

VOBOZBIN=../bin

# build the voronoi diagrams
$VOBOZBIN/voz1b1_serial uniform_grid.raw 0.3 1 UG 2 0 0 0
$VOBOZBIN/voz1b1_serial uniform_grid.raw 0.3 1 UG 2 0 0 1
$VOBOZBIN/voz1b1_serial uniform_grid.raw 0.3 1 UG 2 0 1 0
$VOBOZBIN/voz1b1_serial uniform_grid.raw 0.3 1 UG 2 0 1 1
$VOBOZBIN/voz1b1_serial uniform_grid.raw 0.3 1 UG 2 1 0 0
$VOBOZBIN/voz1b1_serial uniform_grid.raw 0.3 1 UG 2 1 0 1
$VOBOZBIN/voz1b1_serial uniform_grid.raw 0.3 1 UG 2 1 1 0
$VOBOZBIN/voz1b1_serial uniform_grid.raw 0.3 1 UG 2 1 1 1
$VOBOZBIN/voztie 2 UG

# join the zones
$VOBOZBIN/jovoz adjUG.dat volUG.dat zonesUG.dat zonesUG.txt 0.01

# unbind particles
$VOBOZBIN/boz 1.0 0.3 1.0 uniform_grid.raw uniform_grid_vel.raw zonesUG.dat halosUG.dat halosUG.txt 1.2

cat halosUG.txt