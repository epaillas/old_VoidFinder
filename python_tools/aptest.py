import numpy as np
import sys
import os
import glob
import subprocess
from astropy.io import fits
from python_tools.cosmology import Cosmology
from python_tools.galaxycat import GalaxyCatalogue
from python_tools.voidcat import VoidCatalogue
from scipy.spatial import Delaunay
import healpy as hp
from scipy.io import FortranFile
import matplotlib.pyplot as plt

class VoidStatistics:

    def __init__(self, void_file, tracer_file, is_box=True, random_file=None,
                 boss_like=False, pos_cols='0,1,2', box_size=1024.0,
                 omega_m=0.31, h=0.6777, verbose=False, handle=None,
                 ncores=1:

        self.void_file = void_file
        self.tracer_file = tracer_file
        self.random_file = random_file
        self.handle = handle
        self.tracer_unf = self.handle + '.dat.unf'
        self.random_unf = self.handle + '.ran.unf'
        self.void_unf = self.handle + '.voids.unf'
        self.ncores = ncores
        self.is_box = is_box
        self.box_size = box_size
        self.gridmin = -5000
        self.gridmax = 5000
        self.ngrid = 100
        self.dmin = 0
        self.dmax = 3
        self.nbins = 30
        self.min_rvoid = 0
        self.max_rvoid = 200

        # set cosmology
        self.omega_m = omega_m
        self.h = h
        self.cosmo = Cosmology(omega_m=omega_m, h=h)

        pos_cols = [int(i) for i in pos_cols.split(',')]

        self.tracers = GalaxyCatalogue(catalogue_file=tracer_file, is_box=is_box,
            box_size=box_size, randoms=False, boss_like=boss_like, omega_m=omega_m,
            h=h, bin_write=True, output_file=self.tracer_unf, pos_cols=pos_cols,
            zmin=zmin, zmax=zmax)

        self.randoms = GalaxyCatalogue(catalogue_file=random_file, is_box=self.is_box, 
                                       randoms=True, boss_like=boss_like, omega_m=omega_m,
                                       h=h, bin_write=True, output_file=self.random_unf,
                                       pos_cols=pos_cols, zmin=zmin, zmax=zmax)

        self.voids = VoidCatalogue(catalogue_file=void_file, is_box=self.is_box, 
                                   omega_m=omega_m, h=h, bin_write=True,
                                   output_file=self.void_unf, pos_cols=pos_cols)

    def _get_mean_monopole(self, fname, rmin=0, rmax=100):
        
        # read original data
        data= np.genfromtxt(fname, skip_header=1)
        bins = np.genfromtxt(fname, skip_footer=len(data))

        # subsample data by void radius
        data = np.asarray([i for i in data if rmin <= i[0] <= rmax])

        # get mean profile and save data
        mean_profile = np.mean(data[:, 1:], axis=0)
        cout = np.c_[bins, mean_profile]
        fout = self.void_file + '.mean_VG_CCF_monopole'
        fmt = 2*'%10.3f '
        np.savetxt(fout, cout, fmt=fmt)

    
    def VoidGalaxyCCF(self, kind='monopole'):
        if kind == 'monopole':
            self._2PCF_monopole()
        elif kind == 'r-mu':
            self._2PCF_r_mu()
        elif kind == 'sigma-pi':
            self._2PCF_sigma_pi()
        else:
            sys.exit('Correlation kind not recognized. Aborting...')

    def _2PCF_monopole(self):
        '''
        Computes the monopole void-galaxy
        cross-correlation function (bins in
        distance).
        '''
        fout = self.void_file + '.VG_CCF_monopole'
        nbins = 30
        dmin = 0
        dmax = 3

        if self.is_box:
            binpath = sys.path[0] + '/SVF_box/bin/'
            cmd = [binpath + 'vg_ccf_monopole.exe',
                   self.tracer_file,
                   self.void_file,
                   fout,
                   str(self.box_size),
                   str(nbins),
                   str(dmin),
                   str(dmax),
                   ]
        else:
            binpath = sys.path[0] + '/SVF_survey/bin/'
            sys.exit('Not implemented!')

        logfile = self.handle + '_VG_CCF_monopole.log'
        log = open(logfile, "w+")
        subprocess.call(cmd, stdout=log, stderr=log)

        self._get_mean_monopole(fout)

    def _2PCF_r_mu(self):
        '''
        Computes the void-galaxy cross-correlation
        function in bins of r and mu.
        '''

        fout = self.handle + '.VG_CCF_rmu'

        if self.is_box:
            sys.exit('Not implemented!')
        else:
            binpath = sys.path[0] + '/SVF_survey/bin/'
            cmd = [binpath + 'vg_ccf_r_mu.exe',
                   self.tracer_file,
                   self.random_file,
                   self.void_file,
                   fout,
                   str(self.dmin),
                   str(self.dmax),
                   str(self.nbins),
                   str(self.gridmin),
                   str(self.gridmax),
                   str(self.min_rvoid),
                   str(self.max_rvoid)]

        logfile = self.handle + '_vg_ccf_rmu.log'
        log = open(logfile, "w+")
        subprocess.call(cmd, stdout=log, stderr=log)


    def _2PCF_sigma_pi(self):
        '''
        Computes the void-galaxy cross-correlation
        function in bins of sigma and pi.
        '''

        fout = self.handle + '.VG_CCF_spi'

        if self.is_box:
            sys.exit('Not implemented!')
        else:
            binpath = sys.path[0] + '/SVF_survey/bin/'
            cmd = [binpath + 'vg_ccf_sigma_pi.exe',
                   self.tracer_file,
                   self.random_file,
                   self.void_file,
                   fout,
                   str(self.dmin),
                   str(self.dmax),
                   str(self.nbins),
                   str(self.gridmin),
                   str(self.gridmax),
                   str(self.min_rvoid),
                   str(self.max_rvoid)]

        logfile = self.handle + '_vg_ccf_spi.log'
        log = open(logfile, "w+")
        subprocess.call(cmd, stdout=log, stderr=log)


    def bin_centres(self, bin_edges):
        '''
        Find the bin centres
        for an input set of 
        bin edges.
        '''
        return np.asarray([(bin_edges[i] + bin_edges[i + 1])/2
        for i in range(len(bin_edges) - 1)])

    def VoidAbundance(self):

        fout = self.handle + '.void_abundance'

        voids = np.genfromtxt(self.void_file)
        radius = voids[:,3]

        minr = radius.min()
        maxr = radius.max()
        nbins = 10
        bins = np.logspace(np.log10(minr), np.log10(maxr), nbins)
        hist, bin_edges = np.histogram(radius, bins=bins)
       
        if self.is_box:
            norm = (self.box_size**3 * np.diff(np.log10(bin_edges)))
        else:
            norm = 1

        n = hist / norm

        err = np.sqrt(hist) / norm
        errplus = n + err
        errminus = n - err

        bin_centres = self.bin_centres(bin_edges)
        
        cout = np.c_[bin_centres, n, errminus, errplus]
        fmt = 4*'%10.10f '
        np.savetxt(fout, cout)
