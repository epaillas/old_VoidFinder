import click
from python_tools.circularvoids import CircularVoids

@click.command()
@click.option('--tracers', type=str, help='File containing tracers.')
@click.option('--handle', type=str, help='Basename for the output files')
@click.option('--is_box', type=bool, default=True, help='Is the data from a simulation box?')
@click.option('--is_periodic', type=bool, default=True, help='Is the simulation box periodic?')
@click.option('--ncores', type=int, default=1, help='Number of cores to use for parallel tasks.')
@click.option('--steps', type=str, default='1,2,3,4', help='Which steps are to be run. (e.g. 1,2,3).')
@click.option('--pos_cols', type=str, default='0,1,2', help='Indices of columns where tracer positions are stored.')
@click.option('--rvoidmax', type=float, default=50, help='Maximum void radius to search.')
@click.option('--box_size', type=float, default=1024, help='[Periodic box] Size of the simulation box')
@click.option('--delta_voids', type=float, default=0.2, help='Void density threshold.')
@click.option('--randoms', type=str, default='', help='[Survey] File containing randoms.')
@click.option('--mask', type=str, default='', help='[Survey] File containing mask of survey footprint.')
@click.option('--boss_like', type=bool, default=False, help='[Survey] Is the data from BOSS/eBOSS?')
@click.option('--zmin', type=float, default=0.43, help='[Survey] Low redshift cut')
@click.option('--zmax', type=float, default=0.7, help='[Survey] High redshift cut')
@click.option('--omega_m', type=float, default=0.31, help='[Survey] Omega matter.')
@click.option('--h', type=float, default=67.77, help='[Survey] Little h.')
@click.option('--delete_files', type=bool, default=False, help='Delete intermediate files.')

def run_circular_voids(tracers, handle, is_box, ncores, steps, pos_cols,
                        rvoidmax, box_size, randoms, mask, boss_like,
                        zmin, zmax, omega_m, h, is_periodic, delta_voids,
                        delete_files):
    
    voids = CircularVoids(tracer_file=tracers, random_file=randoms, handle=handle,
                           is_box=is_box, box_size=box_size, steps=steps,
                           ncores=ncores, boss_like=boss_like, mask_file=mask,
                           pos_cols=pos_cols, omega_m=omega_m, h=h, delta_voids=delta_voids,
                           rvoidmax=rvoidmax, zmin=zmin, zmax=zmax, is_periodic=is_periodic,
                           delete_files=delete_files)

if __name__ == '__main__':
    run_circular_voids()