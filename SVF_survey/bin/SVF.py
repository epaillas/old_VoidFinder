from galaxycat import GalaxyCatalogue
import numpy as np
import sys
from scipy.spatial import Delaunay



class SphericalVoids:

    def __init__(self, tracer_file, is_box=True, boss_like=False, pos_columns=[0, 1, 2],
                 box_length=1024.0, omega_m=0.31, h=0.6777, mask_file='', z_min=0.43, zmax=0.7,
                 verbose=False, output_folder='', handle=''):

        print('\nRunning the Spherical Void Finder')

        # load tracer information
        self.tracers = GalaxyCatalogue(catalogue_file=tracer_file, is_box=is_box,
        box_length=box_length, randoms=False, boss_like=boss_like, omega_m=omega_m,
        h=h)

    def delaunay_triangulation(self):
        '''
        Make a Delaunay triangulation over
        the cartesian positions of the tracers.
        Returns the vertices of tetrahedra.
        '''
        points = np.hstack([self.tracers.x, self.tracers.y, self.tracers.z])
        triangulation = Delaunay(points)
        simplices = triangulation.simplices.copy()
        vertices = points[simplices]
        return vertices
           


            
            

        

