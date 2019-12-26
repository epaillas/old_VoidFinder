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

class FieldVoids:

    def __init__(self, field_file, centres_file, is_box=True, random_file='',
                 boss_like=False, pos_cols='0,1,2', box_size=1024.0,
                 omega_m=0.31, h=0.6777, mask_file='', zmin=0.43, zmax=0.7,
                 verbose=False, handle='', nside=512, delta_voids=0.2,
                 rvoidmax=50, ncores=1, steps='1,2,3,4', is_periodic=True):

        steps = [int(i) for i in steps.split(',')]
        pos_cols = [int(i) for i in pos_cols.split(',')]

        # file names
        self.handle = handle
        self.field_file = field_file
        self.random_file = random_file
        self.centres_file = centres_file
        self.voids_file = self.handle + '.FVF'
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
        self.cosmo = Cosmology(omega_m=omega_m)


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

            if not self.is_periodic:
                self.remove_edge_voids()

            voids = self.overlap_filter(overlap=0.0)
            voids = self.overlap_filter(overlap=0.2)
            voids = self.overlap_filter(overlap=0.5)


    def concat_files(self, input_files, output_file):
        with open(output_file, 'w+') as outfile:
            for fname in input_files:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)


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
            

    def grow_circles(self, ncores=1):
        '''
        Grow spheres from an input 
        catalogue of centres.
        '''
        print('Proceeding to grow circles...')
        if self.is_box:
            binpath = sys.path[0] + '/SVF_box/bin/'
            self.ngrid = 100
            cmd = ['mpirun', '-np', str(ncores), binpath + 'grow_spheres_2D_DF.exe',
                self.field_file, self.centres_file, self.voids_file, str(self.box_size),
                str(self.delta_voids), str(self.rvoidmax)]
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
        sys.exit('Not implemented!')

    def sort_circles(self, fname='', radius_col=2):
        '''
        Sort an input void catalogue in
        decreasing order of radius.
        '''
        print('Sorting circles by decreasing radius...')

        if fname == '':
            fname = self.recentred_file

        voids = np.genfromtxt(fname)
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


