import click
from python_tools.clusterstatistics import ClusterStatistics

@click.command()
@click.option('--clusters', type=str, required=True, help='File containing clusters.')
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
@click.option('--rcluster_min', type=float, default=0, help='Minimum cluster radius cut')
@click.option('--rcluster_max', type=float, default=300, help='Maximum cluster radius cut')
@click.option('--median_cut', type=bool, default=False, help='Use median cluster radius as a cut.')
@click.option('--dmin', type=float, default=5, help='Minimum radial distance')
@click.option('--dmax', type=float, default=150, help='Maximum radial distance')
def postprocess_clusters(clusters, tracers, randoms, handle, is_box,
                      ncores, box_size, boss_like, pos_cols,
                      velocity, nrbins, rcluster_min, rcluster_max,
                      dmin, dmax, median_cut, is_matter):

    clusters = ClusterStatistics(cluster_file=clusters, tracer_file=tracers, random_file=randoms,
                        handle=handle, is_box=is_box, box_size=box_size, nrbins=nrbins,
                        ncores=ncores, boss_like=boss_like, pos_cols=pos_cols,
                        rcluster_min=rcluster_min, rcluster_max=rcluster_max, dmin=dmin, dmax=dmax)

    if is_matter:
        clusters.ClusterMatterCCF(kind='monopole', median_cut=median_cut)
        clusters.ClusterMatterCCF(kind='r-mu', median_cut=median_cut)

        if velocity:
            clusters.ClusterMatterCCF(kind='los_velocity')
            clusters.ClusterMatterCCF(kind='clustercen_velocity')
    
    else:
        clusters.ClusterGalaxyCCF(kind='monopole', median_cut=median_cut)
        clusters.ClusterGalaxyCCF(kind='r-mu', median_cut=median_cut)
        
        if velocity:
            clusters.ClusterGalaxyCCF(kind='los_velocity')
            clusters.ClusterGalaxyCCF(kind='clustercen_velocity')


if __name__ == '__main__':
    postprocess_clusters()
