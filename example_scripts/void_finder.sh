#!/bin/bash

#-------- Universal settings --------
tracers=$HOME/tracers.dat
handle=$HOME/tracers
is_box=True
ncores=1
steps='1,2,3,4'
pos_cols='0,1,2'
rvoidmax=100

#-------- Periodic box settings --------
box_size=1000

#### Survey settings ####
randoms=$HOME/randoms.dat
mask=
boss_like=False
zmin=0.43
zmax=0.7

#-------- Run (do not modify below this point) --------
python $HOME/code/void_finder/void_finder.py \
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
