import numpy as np
import sys
import os
import glob
import subprocess
from astropy.io import fits
from python_tools.cosmology import Cosmology
from python_tools.galaxycat import ProjectedGalaxyCatalogue
from scipy.spatial import Delaunay
import healpy as hp
from scipy.io import FortranFile
import matplotlib.pyplot as plt

class CircularVoids:

    def __init__(self, tracer_file, is_box=True, random_file='',
                 boss_like=False, pos_cols='0,1,2', box_size=1024.0,
                 omega_m=0.31, h=0.6777, mask_file='', zmin=0.43, zmax=0.7,
                 verbose=False, handle='', nside=512, delta_voids=0.2,
                 rvoidmax=50, ncores=1, steps='1,2,3,4', is_periodic=True,
                 delete_files=False):

        steps = [int(i) for i in steps.split(',')]
        pos_cols = [int(i) for i in pos_cols.split(',')]

        # file names
        self.handle = handle
        self.tracer_file = tracer_file
        self.random_file = random_file
        self.tracer_unf = self.handle + '.dat.unf'
        self.random_unf = self.handle + '.ran.unf'
        self.centres_file = self.handle + '.cen.unf'
        self.voids_file = self.handle + '.SVF'
        self.recentred_file = self.voids_file + '_recen'
        self.mask_file = mask_file
        self.steps = steps
        self.pos_cols = pos_cols
        self.use_guards = True

        # void parameters
        self.delta_voids = delta_voids
        self.rvoidmax = rvoidmax

        # catalog data
        self.is_box = is_box
        self.is_periodic = is_periodic
        self.box_size = box_size
        self.zmin = zmin
        self.zmax = zmax
        self.nside = nside

        # set cosmology
        self.omega_m = omega_m
        self.h = h
        self.cosmo = Cosmology(om_m=omega_m)

        print('handle: ' + self.handle)
        print('tracer_file: ' + self.tracer_file)
        print('centres_file: ' + self.centres_file)
        print('voids_file: ' + self.voids_file)
        print('recentred_file: ' + self.recentred_file)

        if 1 not in steps:
            if not self.is_box:
                if self.mask_file == '':
                    sys.exit('Mask file not provided. Aborting...')
                else:
                    self.mask = hp.read_map(self.mask_file, nest=False, verbose=False)
        else:
        # load tracers and find void centres
            self.tracers = ProjectedGalaxyCatalogue(catalogue_file=tracer_file, is_box=is_box,
            box_size=box_size, randoms=False, boss_like=boss_like, omega_m=omega_m,
            h=h, bin_write=True, output_file=self.tracer_unf, pos_cols=pos_cols,
            zmin=zmin, zmax=zmax)
            
            if self.is_box == False:
                if random_file == '':
                    sys.exit('Random catalogue is missing. Aborting...')
                else:
                    self.randoms = ProjectedGalaxyCatalogue(catalogue_file=random_file, is_box=self.is_box, 
                                                randoms=True, boss_like=boss_like, omega_m=omega_m,
                                                h=h, bin_write=True, output_file=self.random_unf,
                                                pos_cols=pos_cols, zmin=zmin, zmax=zmax)
                    
                if self.mask_file == '':
                    print('No mask file provided. Generating a rough mask...')
                    self.mask = self.make_survey_mask()
                else:
                    self.mask = hp.read_map(self.mask_file, nest=False, verbose=False)
                    
                self.mask_borders = self.get_mask_borders()

            self.get_circumcentres()

        # Grow spheres from the centres found in the previous step
        if 2 in steps:
            voids = self.grow_circles(ncores=ncores)

        # Find better void centres by shifting the original positions
        if 3 in steps:
            voids = self.recentre_circles(ncores=ncores)

        # Sort spheres in decreasing order of radius
        # and remove overlapping spheres
        if 4 in steps:
            voids = self.sort_circles()

            if not self.is_box:
                voids = self.filter_by_volume_fraction(threshold=0.95)
            else:
                if not self.is_periodic:
                    self.remove_edge_voids()

            voids = self.overlap_filter(overlap=0.0)
            voids = self.overlap_filter(overlap=0.2)
            voids = self.overlap_filter(overlap=0.5)

            # save a catalog with sky coordinates
            if not self.is_box:
                self.get_void_skycoords()

        if delete_files:
            os.remove(self.tracer_unf)
            os.remove(self.centres_file)
            os.remove(self.voids_file)
            os.remove(self.recentred_file)

            if not self.is_box:
                os.remove(self.random_unf)
            

    def concat_files(self, input_files, output_file):
        with open(output_file, 'w+') as outfile:
            for fname in input_files:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

    def in_ipynb():
        try:
            cfg = get_ipython().config 
            if cfg['IPKernelApp']['parent_appname'] == 'ipython-notebook':
                return True
            else:
                return False
        except NameError:
            return False

    def remove_edge_voids(self, fname=''):
        print('Removing edge voids...')
        if fname == '':
            fname = self.recentred_file

        voids = np.genfromtxt(fname)
        x = voids[:,0]
        y = voids[:,1]
        radius = voids[:,2]

        condx = np.logical_or(x - radius < 0, x + radius > self.box_size)
        condy = np.logical_or(y - radius < 0, y + radius > self.box_size)

        cond = np.logical_or.reduce((condx, condy))

        voids = voids[~cond]
        np.savetxt(fname, voids)


    def get_periodic_images(self, data):
        '''
        Find the relevant images of a 
        set of points in a box that
        has boundary conditions.
        '''
        images = []
        buffer = self.box_size / 10 # probably an overkill

        for point in data:
            condx = ((self.box_size - buffer) < point[0]) or (point[0] < buffer)
            condy = ((self.box_size - buffer) < point[1]) or (point[1] < buffer)

            if condx and condy:
                shiftx = point[0] + np.copysign(self.box_size, (buffer - point[0]))
                shifty = point[1] + np.copysign(self.box_size, (buffer - point[1]))
                images.append([shiftx, shifty])
            if condx:
                shiftx = point[0] + np.copysign(self.box_size, (buffer - point[0]))
                shifty = point[1]
                images.append([shiftx, shifty])
            if condy:
                shiftx = point[0]
                shifty = point[1] + np.copysign(self.box_size, (buffer - point[1]))
                images.append([shiftx, shifty])

        images = np.asarray(images)
        return images
            

    def make_survey_mask(self):
        sys.exit('Not implemented!')
    
    def get_mask_borders(self):
        sys.exit('Not implemented!')

    def gen_random_sphere(self):
        sys.exit('Not implemented!')

    def gen_guard_particles(self):

        sys.exit('Not implemented!')
     
    def get_circumcentres(self, radius_limit=1000, bin_write=True):
        '''
        Find the centre of the circumspheres
        associated to an input catalogue of
        tetrahedra.
        '''

        print('Finding circumcentres of triangles...')
        cenx, ceny, r = [], [], []
        vertices = self.delaunay_triangulation()
        
        sing = 0
        for tetra in vertices:
            x0, x1, x2 = tetra
            A = []
            B = []
            A.append((x1 - x0).T)
            A.append((x2 - x0).T)
            A = np.asarray(A)
            B = np.sum(A**2, axis=1)
            B = np.asarray(B)
            try:
                C = np.linalg.inv(A).dot(B)
            except:
                sing += 1
                continue
            centre = x0 + 0.5 * C
            radius = 0.5 * np.sqrt(np.sum(C**2))
            if radius < radius_limit:
                cenx.append(centre[0])
                ceny.append(centre[1])
                r.append(radius)

        print('{} singular matrices found.'.format(sing))
                
        cenx = np.asarray(cenx)
        ceny = np.asarray(ceny)

        if not self.is_box:
            sys.exit('Not implemented!')

        else:
            # keep only those centres inside the box
            in_box = ((cenx >= 0) & (cenx <= self.box_size) &
                      (ceny >= 0) & (ceny <= self.box_size)
                     )
            cenx = cenx[in_box]
            ceny = ceny[in_box]

        cenx = cenx.reshape(len(cenx), 1)
        ceny = ceny.reshape(len(ceny), 1)

        cout = np.hstack([cenx, ceny])

        print('{} centres found.'.format(len(cout)))

        if bin_write:
            f = FortranFile(self.centres_file, 'w')
            nrows, ncols = np.shape(cout)
            f.write_record(nrows)
            f.write_record(ncols)
            f.write_record(cout)
            f.close()
        else:
            np.savetxt(self.centres_file, cout)

        self.centres = cout

        return
        

    def delaunay_triangulation(self, guards=False):
        '''
        Make a Delaunay triangulation over
        the cartesian positions of the tracers.
        Returns the vertices of tetrahedra.
        '''
        x = self.tracers.x
        y = self.tracers.y
        points = np.hstack([x, y])

        if self.is_box == False and self.use_guards == True:
            sys.exit('Not implemented!')

        # add periodic images if dealing with a box
        if self.is_box:
            images = self.get_periodic_images(points)
            points = np.vstack([points, images])

        triangulation = Delaunay(points)
        simplices = triangulation.simplices.copy()
        vertices = points[simplices]
        print('{} vertices found.'.format(len(vertices)))
        return vertices

    def grow_circles(self, ncores=1):
        '''
        Grow spheres from an input 
        catalogue of centres.
        '''
        print('Proceeding to grow circles...')
        if self.is_box:
            binpath = sys.path[0] + '/SVF_box/bin/'
            self.ngrid = 100
            cmd = ['mpirun', '-np', str(ncores), binpath + 'grow_spheres_2D.exe',
                self.tracer_unf, self.centres_file, self.voids_file, str(self.box_size),
                str(self.delta_voids), str(self.rvoidmax), str(self.ngrid)]
        else:
            sys.exit('Not implemented!')

        logfile = self.handle + '_grow_circles.log'
        log = open(logfile, "w+")

        subprocess.call(cmd, stdout=log, stderr=log)

        if ncores > 1:
            files = glob.glob(self.voids_file + '.*')
            self.concat_files(input_files=files, output_file=self.voids_file)
            subprocess.call(['rm'] + files)

        voids = np.genfromtxt(self.voids_file)
        return voids

    def recentre_circles(self, ncores=1):
        '''
        Find better centres for an input
        void catalogue.
        '''
        print('Recentring spheres...')

        if self.is_box:
            binpath = sys.path[0] + '/SVF_box/bin/'
            self.ngrid = 100
            cmd = ['mpirun', '-np', str(ncores), binpath + 'recentring_2D.exe',
                self.tracer_unf, self.voids_file, self.recentred_file,
                str(self.box_size), str(self.delta_voids), str(self.rvoidmax),
                str(self.ngrid)]
        else:
            sys.exit('Not implemented!')

        logfile = self.handle + '_recentring.log'
        log = open(logfile, "w+")

        subprocess.call(cmd, stdout=log, stderr=log)

        if ncores > 1:
            files = glob.glob(self.recentred_file + '.*')
            self.concat_files(input_files=files, output_file=self.recentred_file)
            subprocess.call(['rm'] + files)

        voids = np.genfromtxt(self.recentred_file)
        return voids

    def sort_circles(self, fname='', radius_col=2):
        '''
        Sort an input void catalogue in
        decreasing order of radius.
        '''
        print('Sorting circles by decreasing radius...')

        if fname == '':
            fname = self.recentred_file

        voids = np.genfromtxt(fname)
        if len(voids) < 1:
            voids = voids[np.argsort(voids[:, radius_col])]
        voids = voids[::-1]
        
        fmt = 3*'%10.3f ' + '%10i ' + '%10.3f '
        np.savetxt(fname, voids, fmt=fmt)

        return voids

    def overlap_filter(self, overlap=0.0):

        self.filtered_file = self.recentred_file + '_ovl{}'.format(str(overlap))

        if self.is_box:
            binpath = sys.path[0] + '/SVF_box/bin/'
            self.ngrid = 100
            cmd = [binpath + 'overlapping_2D.exe', self.recentred_file, self.filtered_file,
                   str(self.box_size), str(overlap), str(self.ngrid)]
        else:
            sys.exit('Not implemented!')

        logfile = self.handle + '_overlapping.log'
        log = open(logfile, "w+")

        subprocess.call(cmd, stdout=log, stderr=log)

        voids = np.genfromtxt(self.filtered_file)
        return voids


    def filter_by_volume_fraction(self, threshold=0.95):
        sys.exit('Not implemented!')


    def get_void_volume_fraction(self, fname='', pos_cols=[0,1,2],
                                 radius_col=3):
        sys.exit('Not implemented!')

    def get_void_skycoords(self, fname=''):
        sys.exit('Not implemented!')