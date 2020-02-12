#!/bin/bash -f

# this script was created by vozinit, and modified to 
# correct the path and do all the work
# Rick Wagner - 15JAN08

VOBOZBIN=../bin

# build the voronoi diagrams
$VOBOZBIN/voz1b1_serial single_halo.raw 0.3 1 SH 2 0 0 0
$VOBOZBIN/voz1b1_serial single_halo.raw 0.3 1 SH 2 0 0 1
$VOBOZBIN/voz1b1_serial single_halo.raw 0.3 1 SH 2 0 1 0
$VOBOZBIN/voz1b1_serial single_halo.raw 0.3 1 SH 2 0 1 1
$VOBOZBIN/voz1b1_serial single_halo.raw 0.3 1 SH 2 1 0 0
$VOBOZBIN/voz1b1_serial single_halo.raw 0.3 1 SH 2 1 0 1
$VOBOZBIN/voz1b1_serial single_halo.raw 0.3 1 SH 2 1 1 0
$VOBOZBIN/voz1b1_serial single_halo.raw 0.3 1 SH 2 1 1 1
$VOBOZBIN/voztie 2 SH

# join the zones
$VOBOZBIN/jovoz adjSH.dat volSH.dat zonesSH.dat zonesSH.txt 0.01

# unbind particles
$VOBOZBIN/boz 1.0 0.3 1.0 single_halo.raw single_halo_vel.raw zonesSH.dat halosSH.dat halosSH.txt 1.2

cat halosSH.txt