import numpy as np
import sys
import os
import glob
import subprocess
from astropy.io import fits
from python_tools.cosmology import Cosmology
from python_tools.galaxycat import GalaxyCatalogue
from scipy.spatial import Delaunay
import healpy as hp
from scipy.io import FortranFile
import matplotlib.pyplot as plt

class VoidStatistics:

    def __init__(self, void_file, tracer_file, is_box=True, random_file=None,
                 boss_like=False, pos_cols='0,1,2', box_size=1024.0,
                 omega_m=0.31, h=0.6777, verbose=False, handle=None,
                 ncores=1, nrbins=60, rvoid_min=0, rvoid_max=500,
                 dmin=0, dmax=3, gridmin=-5000, gridmax=5000, ngrid=100):

        self.void_file = void_file
        self.tracer_file = tracer_file
        self.random_file = random_file
        self.handle = handle
        self.ncores = ncores
        self.is_box = is_box
        self.box_size = box_size
        self.gridmin = gridmin
        self.gridmax = gridmax
        self.ngrid = ngrid
        self.dmin = dmin
        self.dmax = dmax
        self.nrbins = nrbins
        self.rvoid_min = rvoid_min
        self.rvoid_max = rvoid_max
        self.is_matter = False

        # set cosmology
        self.omega_m = omega_m
        self.h = h
        self.cosmo = Cosmology(om_m=omega_m, h=h)

        pos_cols = [int(i) for i in pos_cols.split(',')]

   
    def VoidMatterCCF(self, kind='monopole', median_cut=False):
        self.is_matter = True

        if median_cut:
            data = np.genfromtxt(self.void_file)
            rv = data[:,3]
            self.rvoid_min = np.median(rv)
            self.rvoid_max = 500

        if kind == 'monopole':
            self._2PCF_monopole()
        elif kind == 'r-mu':
            self._2PCF_r_mu()
        elif kind == 'sigma-pi':
            self._2PCF_sigma_pi()
        elif kind == 'los_velocity':
            self._2PCF_los_velocity()
        elif kind == 'voidcen_velocity':
            self._2PCF_voidcen_velocity()
        else:
            sys.exit('Correlation kind not recognized. Aborting...')


    def VoidGalaxyCCF(self, kind='monopole', median_cut=False):

        if median_cut:
            data = np.genfromtxt(self.void_file)
            rv = data[:,3]
            self.rvoid_min = np.median(rv)
            self.rvoid_max = 500

        if kind == 'monopole':
            self._2PCF_monopole()
        elif kind == 'r-mu':
            self._2PCF_r_mu()
        elif kind == 'sigma-pi':
            self._2PCF_sigma_pi()
        elif kind == 'los_velocity':
            self._2PCF_los_velocity()
        elif kind == 'voidcen_velocity':
            self._2PCF_voidcen_velocity()
        else:
            sys.exit('Correlation kind not recognized. Aborting...')

    def _2PCF_monopole(self):
        '''
        Computes the monopole void-galaxy
        cross-correlation function (bins in
        distance).
        '''
        if self.is_matter: 
            fout = self.handle + '.VM_CCF_monopole'
            logfile = self.handle + '_vm_ccf_monopole.log'
        else:
            fout = self.handle + '.VG_CCF_monopole'
            logfile = self.handle + '_vg_ccf_monopole.log'

        if self.is_box:
            binpath = sys.path[0] + '/SVF_box/bin/'
            cmd = [binpath + 'vg_ccf_monopole.exe',
                   self.tracer_file,
                   self.void_file,
                   fout,
                   str(self.box_size),
                   str(self.dmin),
                   str(self.dmax),
                   str(self.nrbins),
                   str(self.rvoid_min),
                   str(self.rvoid_max),
                   str(self.ngrid)]
        else:
            binpath = sys.path[0] + '/SVF_survey/bin/'
            sys.exit('Not implemented!')

        log = open(logfile, "w+")
        subprocess.call(cmd, stdout=log, stderr=log)


    def _2PCF_r_mu(self):
        '''
        Computes the void-galaxy cross-correlation
        function in bins of r and mu.
        '''

        if self.is_matter: 
            fout = self.handle + '.VM_CCF_rmu'
            logfile = self.handle + '_vm_ccf_rmu.log'

        else:
            logfile = self.handle + '_vg_ccf_rmu.log'
            fout = self.handle + '.VG_CCF_rmu'

        if self.is_box:
            binpath = sys.path[0] + '/SVF_box/bin/'
            cmd = [binpath + 'vg_ccf_r_mu.exe',
                   self.tracer_file,
                   self.void_file,
                   fout,
                   str(self.box_size),
                   str(self.dmin),
                   str(self.dmax),
                   str(self.nrbins),
                   str(self.rvoid_min),
                   str(self.rvoid_max),
                   str(self.ngrid)]
        else:
            binpath = sys.path[0] + '/SVF_survey/bin/'
            cmd = [binpath + 'vg_ccf_r_mu.exe',
                   self.tracer_file,
                   self.random_file,
                   self.void_file,
                   fout,
                   str(self.dmin),
                   str(self.dmax),
                   str(self.nrbins),
                   str(self.gridmin),
                   str(self.gridmax),
                   str(self.rvoid_min),
                   str(self.rvoid_max)]

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
                   str(self.nrbins),
                   str(self.gridmin),
                   str(self.gridmax),
                   str(self.rvoid_min),
                   str(self.rvoid_max)]

        logfile = self.handle + '_vg_ccf_spi.log'
        log = open(logfile, "w+")
        subprocess.call(cmd, stdout=log, stderr=log)

    def _2PCF_los_velocity(self):
        '''
        Computes the void-galaxy cross-correlation
        function in bins of r and mu.
        '''

        if self.is_matter: 
            fout = self.handle + '.VM_CCF_losvel'
            logfile = self.handle + '_vm_ccf_losvel.log'
        else:
            fout = self.handle + '.VG_CCF_losvel'
            logfile = self.handle + '_vg_ccf_losvel.log'

        if self.is_box:
            binpath = sys.path[0] + '/SVF_box/bin/'
            cmd = [binpath + 'vg_ccf_los_velocity.exe',
                   self.tracer_file,
                   self.void_file,
                   fout,
                   str(self.box_size),
                   str(self.dmin),
                   str(self.dmax),
                   str(self.nrbins),
                   str(self.rvoid_min),
                   str(self.rvoid_max),
                   str(self.ngrid)]
        else:
            sys.exit('Not implemented...')

        log = open(logfile, "w+")
        subprocess.call(cmd, stdout=log, stderr=log)


    def _2PCF_voidcen_velocity(self):
        '''
        Computes the void-galaxy cross-correlation
        function in bins of r and mu.
        '''

        if self.is_matter: 
            fout = self.handle + '.VM_CCF_voidcenvel'
            logfile = self.handle + '_vm_ccf_voidcenvel.log'
        else:
            fout = self.handle + '.VG_CCF_voidcenvel'
            logfile = self.handle + '_vg_ccf_voidcenvel.log'

        if self.is_box:
            binpath = sys.path[0] + '/SVF_box/bin/'
            cmd = [binpath + 'vg_ccf_voidcen_velocity.exe',
                   self.tracer_file,
                   self.void_file,
                   fout,
                   str(self.box_size),
                   str(self.dmin),
                   str(self.dmax),
                   str(self.nrbins),
                   str(self.rvoid_min),
                   str(self.rvoid_max),
                   str(self.ngrid)]
        else:
            sys.exit('Not implemented...')

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
        nrbins = 10
        bins = np.logspace(np.log10(minr), np.log10(maxr), nrbins)
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
