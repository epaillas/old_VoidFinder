#!/bin/bash -f

# this script was created by vozinit, and modified to 
# correct the path and do all the work
# Rick Wagner - 15JAN08

VOBOZBIN=../bin

# build the voronoi diagrams
$VOBOZBIN/voz1b1_serial single_halo_boundary.raw 0.3 1 SHB 2 0 0 0
$VOBOZBIN/voz1b1_serial single_halo_boundary.raw 0.3 1 SHB 2 0 0 1
$VOBOZBIN/voz1b1_serial single_halo_boundary.raw 0.3 1 SHB 2 0 1 0
$VOBOZBIN/voz1b1_serial single_halo_boundary.raw 0.3 1 SHB 2 0 1 1
$VOBOZBIN/voz1b1_serial single_halo_boundary.raw 0.3 1 SHB 2 1 0 0
$VOBOZBIN/voz1b1_serial single_halo_boundary.raw 0.3 1 SHB 2 1 0 1
$VOBOZBIN/voz1b1_serial single_halo_boundary.raw 0.3 1 SHB 2 1 1 0
$VOBOZBIN/voz1b1_serial single_halo_boundary.raw 0.3 1 SHB 2 1 1 1
$VOBOZBIN/voztie 2 SHB

# join the zones
$VOBOZBIN/jovoz adjSHB.dat volSHB.dat zonesSHB.dat zonesSHB.txt 0.01

# unbind particles
$VOBOZBIN/boz 1.0 0.3 1.0 single_halo_boundary.raw single_halo_boundary_vel.raw zonesSHB.dat halosSHB.dat halosSHB.txt 1.2

cat halosSHB.txt