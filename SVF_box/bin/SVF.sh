#!/bin/bash
set -e

# Python path
python="/anaconda3/bin/python"

bin_dir="/Users/epaillas/data/new_svf/SVF_box_delaunay/bin"
tracers="HOD_GR_B1024_Box1_z0.0.dat"
out_basename="HOD_GR_B1024_Box1_z0.0"
boxsize=1024
buffer=20
delta=0.2
ngrid=64
rvoidmax=100
ncores=80

# Delaunay triangulation
# 1) input_tracers
# 2) output_vertices
# 3) boxsize
# 4) buffer_size
$python $bin_dir/delaunay_triangulation.py $tracers $out_basename".ver" $boxsize $buffer

# Circumcentres
# 1) input_vertices
# 2) output_centres
# 3) boxsize
$bin_dir/circumcentre.exe $out_basename".ver" $out_basename".cen" $boxsize

# Grow spheres
# 1) input_tracers
# 2) input_centres
# 3) output_voids
# 4) boxsize
# 5) density_threshold
# 6) rvoidmax
# 7) ngrid
mpirun -np $ncores $bin_dir/grow_spheres.exe $tracers $out_basename".cen" $out_basename".SVF" $boxsize $delta $rvoidmax $ngrid

if [ $ncores -gt 1 ]; then
  cat $out_basename".SVF.*" > $out_basename".SVF"
  rm $out_basename".SVF.*"
fi

# Recentring
# 1) input_tracers
# 2) input_centres
# 3) output_voids
# 4) boxsize
# 5) density_threshold
# 6) rvoidmax
# 7) ngrid
mpirun -np $ncores $bin_dir/recentring.exe $tracers $out_basename".SVF" $out_basename".SVF.recen" $boxsize $delta $rvoidmax $ngrid

if [ $ncores -gt 1 ]; then
  cat $out_basename".SVF.recen.*" > $out_basename".SVF.recen"
  rm $out_basename".SVF.recen.*"
fi

# Sorting
# 1) input_voids
$python $bin_dir/sort.py $out_basename".SVF.recen"

# Overlapping
# 1) input_voids
# 2) output_voids
# 3) boxsize
# 4) overlap
# 5) ngrid
$bin_dir/overlapping.exe $out_basename".SVF.recen" $out_basename".SVF.recen.ovl0.0" $boxsize 0.0 $ngrid
$bin_dir/overlapping.exe $out_basename".SVF.recen" $out_basename".SVF.recen.ovl0.2" $boxsize 0.2 $ngrid
$bin_dir/overlapping.exe $out_basename".SVF.recen" $out_basename".SVF.recen.ovl0.5" $boxsize 0.5 $ngrid
