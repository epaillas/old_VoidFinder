#!/bin/bash -f

# Mark Neyrinck, March 21 2009

QHULLBIN=../qhull
VOBOZBIN=../bin

# 1024 uniform Poisson-distributed points in a unit box centered at the origin,
# surrounded by 1024 buffer points in a shell of thickness 0.1

$QHULLBIN/rbox 1024 D3 B0.5 | tail -1024 > inner.pos
$QHULLBIN/rbox 1024 D3 B0.6 W0.1 | tail -1024 > buffer.pos
echo 2048 1024 > header.txt
cat header.txt inner.pos buffer.pos > PVBUF.pos
rm header.txt inner.pos buffer.pos

echo
echo ****************
echo * Tessellating *
echo ****************

$VOBOZBIN/vozisol PVBUF 1. n
echo *********************
echo * Void zone-finding *
echo *********************
# change 'v' to 'c' for cluster-finding
$VOBOZBIN/jozov v PVBUF 0 0

echo ****************************************
echo * Python routines to modify void lists *
echo ****************************************
export PYTHONPATH=.:../bin
python poisson_vloidsters.py
