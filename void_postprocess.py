import click
from python_tools.voidstatistics import VoidStatistics

@click.command()
@click.option('--voids', type=str, required=True, help='File containing voids.')
@click.option('--tracers', type=str, help='File containing tracers.')
@click.option('--randoms', type=str, help='File containing randoms.')
@click.option('--handle', type=str, required=True, help='Basename for the output files')
@click.option('--is_box', type=bool, default=True, help='Is the data from a simulation box?')
@click.option('--is_matter', type=bool, default=False, help='Are the tracers DM particles?')
@click.option('--boss_like', type=bool, default=False, help='Is the data from BOSS/eBOSS?')
@click.option('--ncores', type=int, default=1, help='Number of cores to use for parallel tasks.')
@click.option('--box_size', type=float, default=1024, help='Size of the simulation box (only used if is_box is True)')
@click.option('--pos_cols', type=str, default='1,2,3', help='Indices of columns where tracer positions are stored.')
@click.option('--velocity', type=bool, default=False, help='Include velocity statistics?')
@click.option('--nrbins', type=int, default=60, help='Number of radial bins')
@click.option('--rvoid_min', type=float, default=0, help='Minimum void radius cut')
@click.option('--rvoid_max', type=float, default=300, help='Maximum void radius cut')
@click.option('--median_cut', type=bool, default=False, help='Use median void radius as a cut.')
@click.option('--dmin', type=float, default=5, help='Minimum radial distance')
@click.option('--dmax', type=float, default=150, help='Maximum radial distance')
def postprocess_voids(voids, tracers, randoms, handle, is_box,
                      ncores, box_size, boss_like, pos_cols,
                      velocity, nrbins, rvoid_min, rvoid_max,
                      dmin, dmax, median_cut, is_matter):

    voids = VoidStatistics(void_file=voids, tracer_file=tracers, random_file=randoms,
                        handle=handle, is_box=is_box, box_size=box_size, nrbins=nrbins,
                        ncores=ncores, boss_like=boss_like, pos_cols=pos_cols,
                        rvoid_min=rvoid_min, rvoid_max=rvoid_max, dmin=dmin, dmax=dmax)

    if is_matter:
        voids.VoidMatterCCF(kind='monopole', median_cut=median_cut)
        voids.VoidMatterCCF(kind='r-mu', median_cut=median_cut)

        if velocity:
            voids.VoidMatterCCF(kind='los_velocity')
            voids.VoidMatterCCF(kind='voidcen_velocity')
    
    else:
        voids.VoidGalaxyCCF(kind='monopole', median_cut=median_cut)
        voids.VoidGalaxyCCF(kind='r-mu', median_cut=median_cut)
        
        if velocity:
            voids.VoidGalaxyCCF(kind='los_velocity')
            voids.VoidGalaxyCCF(kind='voidcen_velocity')


if __name__ == '__main__':
    postprocess_voids()
