#!/bin/bash

#-------- Universal settings --------
voids=$HOME/voids.dat
tracers=$HOME/tracers.unf
handle=$HOME/tracers
is_box=True
ncores=1
pos_cols='0,1,2'

#-------- Periodic box settings --------
box_size=1000

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
