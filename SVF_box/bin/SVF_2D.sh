#!/bin/bash
set -e

# Python path
python="/anaconda3/bin/python"

bin_dir="/Users/epaillas/data/new_svf/SVF_box_delaunay/bin"
tracers="HOD_GR_B1024_Box1_z0.0.dat"
out_basename="HOD_GR_B1024_Box1_z0.0"
boxsize=1024
buffer=20
delta=0.4
ngrid=64
rvoidmax=20

# Delaunay triangulation
# 1) input_tracers
# 2) output_vertices
# 3) boxsize
# 4) buffer_size
$python $bin_dir/delaunay_triangulation_2D.py $tracers $out_basename".ver_2D" $boxsize $buffer

# Circumcentres
# 1) input_vertices
# 2) output_centres
# 3) boxsize
$bin_dir/circumcentre_2D.exe $out_basename".ver_2D" $out_basename".cen_2D" $boxsize

# Grow spheres
# 1) input_tracers
# 2) input_centres
# 3) output_voids
# 4) boxsize
# 5) density_threshold
# 6) rvoidmax
# 7) ngrid
$bin_dir/grow_spheres_2D.exe $tracers $out_basename".cen_2D" $out_basename".SVF_2D" $boxsize $delta $rvoidmax $ngrid

if [ $ncores -gt 1 ]; then
  cat $out_basename".SVF_2D.*" > $out_basename".SVF_2D"
  rm $out_basename".SVF_2D.*"
fi

# Recentring
# 1) input_tracers
# 2) input_centres
# 3) output_voids
# 4) boxsize
# 5) density_threshold
# 6) rvoidmax
# 7) ngrid
$bin_dir/recentring_2D.exe $tracers $out_basename".SVF_2D" $out_basename".SVF_2D.recen" $boxsize $delta $rvoidmax $ngrid

if [ $ncores -gt 1 ]; then
  cat $out_basename".SVF_2D.recen.*" > $out_basename".SVF_2D.recen"
  rm $out_basename".SVF_2D.recen.*"
fi

# Sorting
# 1) input_voids
$python $bin_dir/sort_2D.py $out_basename".SVF_2D.recen"

# Overlapping
# 1) input_voids
# 2) output_voids
# 3) boxsize
# 4) overlap
# 5) ngrid
$bin_dir/overlapping_2D.exe $out_basename".SVF_2D.recen" $out_basename".SVF_2D.recen.ovl0.0" $boxsize 0.0 $ngrid
$bin_dir/overlapping_2D.exe $out_basename".SVF_2D.recen" $out_basename".SVF_2D.recen.ovl0.2" $boxsize 0.2 $ngrid
$bin_dir/overlapping_2D.exe $out_basename".SVF_2D.recen" $out_basename".SVF_2D.recen.ovl0.5" $boxsize 0.5 $ngrid
