#!/bin/bash

#-------- Universal settings --------
tracers=$HOME/tracers.dat
handle=$HOME/tracers
is_box=True

#-------- Periodic box settings --------
box_size=1000


#-------- Run (do not modify below this point) --------
python $HOME/code/void_finder/zobov_void_finder.py \
--tracers "$tracers" \
--handle "$handle" \
--is_box "$is_box" \
--box_size "$box_size"
