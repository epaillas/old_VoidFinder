#!/bin/bash

#-------- Universal settings --------
voids=$HOME/voids.dat
tracers=$HOME/tracers.unf
handle=$HOME/tracers
is_box=True
ncores=1
pos_cols='0,1,2'
nrbins=30
rvoid_min=30
rvoid_max=300
dmin=5
dmax=150
median_cut=False

#-------- Periodic box settings --------
box_size=1500
velocity=True

#### Survey settings ####
randoms=$HOME/randoms.unf
boss_like=True

#-------- Run (do not modify below this point) --------
python $HOME/code/void_finder/void_postprocess.py \
--voids "$voids" \
--tracers "$tracers" \
--handle "$handle" \
--is_box "$is_box" \
--ncores "$ncores" \
--pos_cols "$pos_cols" \
--box_size "$box_size" \
--randoms "$randoms" \
--boss_like "$boss_like" \
--velocity "$velocity" \
--nrbins "$nrbins" \
--rvoid_min "$rvoid_min" \
--rvoid_max "$rvoid_max" \
--dmin "$dmin" \
--dmax "$dmax" \
--median_cut "$median_cut"

