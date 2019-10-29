#!/bin/bash

#-------- Universal settings --------
voids=$HOME/voids.dat
tracers=$HOME/tracers.unf
handle=$HOME/tracers
is_box=True
ncores=1
steps='1,2,3,4'
pos_cols='0,1,2'
rvoidmax=100

#-------- Periodic box settings --------
box_size=1000

#### Survey settings ####
randoms=$HOME/randoms.unf
mask=
boss_like=True
zmin=0.43
zmax=0.7

#-------- Run (do not modify below this point) --------
python $HOME/code/void_finder/void_postprocess.py \
--voids "$voids" \
--tracers "$tracers" \
--handle "$handle" \
--is_box "$is_box" \
--ncores "$ncores" \
--steps "$steps" \
--pos_cols "$pos_cols" \
--rvoidmax "$rvoidmax" \
--box_size "$box_size" \
--randoms "$randoms" \
--mask "$mask" \
--boss_like "$boss_like" \
--zmin "$zmin" \
--zmax "$zmax" \
