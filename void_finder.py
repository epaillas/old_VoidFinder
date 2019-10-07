import click
from python_tools.sphericalvoids import SphericalVoids

@click.command()
@click.option('--tracers', type=str, required=True, help='File containing tracers.')
@click.option('--randoms', type=str, help='File containing randoms.')
@click.option('--handle', type=str, required=True, help='Basename for the output files')
@click.option('--is_box', type=bool, default=True, help='Is the data from a simulation box? (True or False')
@click.option('--boss_like', type=bool, default=False, help='Is the data from BOSS/eBOSS?')
@click.option('--ncores', type=int, default=1, help='Number of cores to use for parallel tasks.')
@click.option('--box_size', type=float, default=1024, help='Size of the simulation box (only used if is_box is True)')
@click.option('--steps', type=str, default='1,2,3,4', help='Which steps are to be run. (e.g. 1,2,3).')
@click.option('--mask_file', type=str, help='File containing HEALPix mask of survey footprint.')
@click.option('--pos_cols', type=str, default='1,2,3', help='Indices of columns where tracer positions are stored.')
def find_voids(tracers, randoms, handle, is_box, ncores,
                        box_size, steps, boss_like, mask_file,
                        pos_cols):

    voids = SphericalVoids(tracer_file=tracers, random_file=randoms, handle=handle,
                        is_box=is_box, box_size=box_size, steps=steps, ncores=ncores,
                        boss_like=boss_like, mask_file=mask_file, pos_cols=pos_cols)


if __name__ == '__main__':
    find_voids()