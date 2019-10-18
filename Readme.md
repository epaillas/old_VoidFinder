### Void Finder ###

Allows you to indentify spherical voids in a large-scale distribution
of galaxies, haloes or particles. It works for periodic cosmological boxes and
observational surveys, such as SDSS BOSS. In future realeases, multiple void 
finder algorithms will be implemented. The core of the algorithm is written
in Fortran 90, and everything is wrapped under a Python API.

To compile the files, type *make* under each void finder's directory. For example:

```
cd SVF_box
make
cd ../SVF_survey
make
```


To set everything up for void identification, copy and edit *example_script.sh*
under your working directory. Here are the different options that need to be 
tuned before running the code (which can also by accessed by running python
void_finder.py --help):

Options:
* tracers TEXT       File containing tracers.
* handle TEXT        Basename for the output files
* is_box BOOLEAN     Is the data from a simulation box?
* ncores INTEGER     Number of cores to use for parallel tasks.
* steps TEXT         Which steps are to be run. (e.g. 1,2,3).
* pos_cols TEXT      Indices of columns where tracer positions are stored.
* rvoidmax FLOAT     Maximum void radius to search.
* box_size FLOAT     [Periodic box] Size of the simulation box
* randoms TEXT       [Survey] File containing randoms.
* mask TEXT          [Survey] File containing mask of survey footprint.
* boss_like BOOLEAN  [Survey] Is the data from BOSS/eBOSS?
* zmin FLOAT         [Survey] Low redshift cut
* zmax FLOAT         [Survey] High redshift cut
* help               Show this message and exit.

Upon configuring the file, simply run it as

```
sh example_script.sh
```
