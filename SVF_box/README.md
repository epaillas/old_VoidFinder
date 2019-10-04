# Spherical Void Finder (periodic box version)

This void finding algorithm searches for underdense spherical regions in a periodic cosmological simulation volume. The code
is written in a mixture of Fortran 90 and Python 3, and can be run in parallel using MPI. It can be run on the dark matter
density field, as well as on a biased distribution of tracers, such as haloes or galaxies.

## Pipeline summary

1) A Delaunay triangulation is computed from tracer positions, and the circumcentres associated to each Delaunay
tetrahedron are identified. 
    
2) Starting from a fixed radius R_max around each centre, spheres are shrank until their integrated tracer number
density is equal or lower than a specified fraction of the mean density (commonly people use 20 per cent of the mean).
The radius at which this condition is satisfied is stored as the void radius. The value of R_max is chosen
such that it overestimates the expected radius of the largest void in the sample.
    
3) Given that some of the void centres from the previous step might be shifted with respect to the local density minimum,
the original centres are displaced in a random-walk fashion. In each random step, procedure 2) is repeated, updating the
centre of the void only if the new radius is larger than the previous one. In total, 128 random steps are tried for each
void, and the exploration is limited to regions within the radius of the original void centre from 2).
   
4) We rank the voids by decreasing order of radius, and we check for overlapping voids. If two voids have centres that lie
closer than 50 per cent of the sum of their radii, only the largest one is kept in the catalogue. Allowing a degree of
overlap between the voids gives some freedom to the algorithm to trace non-spherical underdensities in the cosmic web and
increases the statistics, although too much overlap may result in unwanted covariances in the results.

## Compilation

Under the main directory, simply type `make`. This will compile the files under the `src/` directory and copy the executables to the `bin/` directory. If you find problems, you might need to change your Fortran and MPI compilers in the Makefile at `src/`.
## Running the code

The pipeline is divided in individual codes that take care of each of the steps outlined in the summary above. If you want to run the full pipeline, the codes should be executed as follows (required arguments are specified for each case):

1. **delaunay_triangulation.py**: Performs the Delaunay triangulation on the tracers.

    1) input_tracers: Name of the input tracers file. It should be an ascii file with first 3 columns having the tracer positions.
    2) output_vertices: Name of the output file that stores vertices of the triangulation.
    3) boxsize: Size of the simulation box

2. **circumcentre.exe**: Find the circumcentres associated to each Delaunay tetrahedron.

    1) input_vertices: Name of the input vertices file obtained in 1).
    2) output_centres: Name of the output file that stores the circumcentres of each tetrahedron.
    3) boxsize: Size of the simulation box.

3. **grow_spheres.exe**: Finds voids around each Delaunay circumcentres.
    1) input_tracers: Name of the input tracers file (same format as in 1)).
    2) inout_centres: Name of the centres file from 2).
    3) output_voids: Name of the output file that stores the void catalogue.
    4) boxsize: Size of the simulation box
    5) density_threshold: Integrated density threshold that defines the void (e.g. 0.2 to use 20% of the mean).
    6) rvoidmax: The maximum radius from which to shrink the spheres. It should be an over-estimation of the expected
    largest void radius in the simulation.
    7) ngrid: Size of the regular grid that is going to be used for the calculations. It should be a divisor of the
    size of the box (e.g. if the box is 1024 Mpc on a side, ngrid could be 64, 128 etc). Having an ngrid that divides 
    the box in roughly 10 parts is recommended.

4. **recentring.exe**: Re-recentres the voids to maximize their correspondence to local density minima.
    1) input_tracers: Name of the input tracers file (same format as in 1)).
    2) inout_centres: Name of the void catalogue from 3).
    3) output_voids: Name of the output file that stores the re.centr-d void catalogue.
    4) boxsize: Size of the simulation box
    5) density_threshold: Integrated density threshold that defines the void (e.g. 0.2 to use 20% of the mean).
    6) rvoidmax: The maximum radius from which to shrink the spheres.
    7) ngrid: Size of the regular grid that is going to be used for the calculations.

5. **sort.py**: Sort voids in decreasing order of radius.
    1) input_voids: name of the void catalogue obtained in 4).

6. **overlapping.exe**: Filters the voids according to how much they overlap with each other.
    1) input_voids: Name of the sorted void catalogue obtained in 5).
    2) output_voids: Name of the output file that stores the filtered void catalogue.
    3) boxsize: Size of the simulation box.
    4) overlap: Percentage of overlap that is allowed for the voids (e.g. 0.5 for 50% overlap).
    5) ngrid: Size of the regular grid that is going to be used for the calculations.
    
    
## Running in parallel

Both **grow_spheres.exe** and **recentring.exe** can be run in parallel using MPI. To do so, run the codes as `mpirun -np $ncores grow_spheres.exe (or recentring.exe)`, where `$ncores` is the number of cores that you want to use.
