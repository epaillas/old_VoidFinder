import click
from python_tools.zobovvoids import ZobovVoids

@click.command()
@click.option('--tracers', type=str, help='File containing tracers.')
@click.option('--handle', type=str, help='Basename for the output files')
@click.option('--is_box', type=bool, default=True, help='Is the data from a simulation box?')
@click.option('--box_size', type=float, default=1024, help='[Periodic box] Size of the simulation box')

def run_zobov_voids(tracers, handle, is_box,
                    box_size):
    
    voids = ZobovVoids(tracer_file=tracers, handle=handle,
                           is_box=is_box, box_size=box_size)

if __name__ == '__main__':
    run_zobov_voids()