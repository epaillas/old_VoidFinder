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
                 boss_like=False, pos_cols='1,2,3', box_size=1024.0,
                 omega_m=0.31, h=0.6777, verbose=False, handle=None,
                 ncores=1):

        self.void_file = void_file
        self.tracer_file = tracer_file
        self.random_file = random_file
        self.handle = handle
        self.ncores = ncores
        self.is_box = is_box
        self.box_size = box_size
        self.gridmin = -5000
        self.gridmax = 5000
        self.ngrid = 100

        # set cosmology
        self.omega_m = omega_m
        self.h = h
        self.cosmo = Cosmology(omega_m=omega_m)

        pos_cols = [int(i) for i in pos_cols.split(',')]

    def VoidGalaxyCCF(self, kind='monopole'):
        if kind == 'monopole':
            self._2PCF_monopole()
        elif kind == 'r-mu':
            self._2PCF_r_mu()
        elif kind == 'sigma-pi'():
            self._2PCF_sigma_pi()
        else:
            sys.exit('Correlation kind not recognized. Aborting...')

    def _2PCF_monopole(self):
        '''
        Computes the monopole void-galaxy
        cross-correlation function (bins in
        distance).
        '''
        fout = self.handle + '_vg_ccf_monopole.dat'
        nbins = 30
        dmin = 0
        dmax = 3

        if self.is_box:
            binpath = sys.path[0] + '/SVF_box/bin/'
            cmd = [binpath + 'vg_ccf_monopole',
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

        logfile = self.handle + '_vg_ccf_monopole.log'
        log = open(logfile, "w+")
        subprocess.call(cmd, stdout=log, stderr=log)

    def _2PCF_r_mu(self):
        '''
        Computes the void-galaxy cross-correlation
        function in bins of r and mu.
        '''

        fout = handle + '_vg_ccf_rmu.dat'

        if self.is_box:
            sys.exit('Not implemented!')
        else:
            binpath = sys.path[0] + '/SVF_survey/bin/'
            cmd = [binpath + 'vg_ccf_r_mu',
                   self.tracer_unf,
                   self.random_unf,
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

        fout = handle + '_vg_ccf_rmu.dat'

        if self.is_box:
            sys.exit('Not implemented!')
        else:
            binpath = sys.path[0] + '/SVF_survey/bin/'
            cmd = [binpath + 'vg_ccf_sigma_pi.exe',
                   self.tracer_unf,
                   self.random_unf,
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


    def bin_centres(self, bin_edges):
        '''
        Find the bin centres
        for an input set of 
        bin edges.
        '''
        return np.asarray([(bin_edges[i] + bin_edges[i + 1])/2
        for i in range(len(bin_edges) - 1)])

    def plot_void_abundance(self, nbins=15, error=True, linestyle='-',
                            linewidth=2.0, color='k'):
        '''
        Plot the distribution 
        of void sizes.
        '''

        voids = np.genfromtxt(self.void_file)
        radius = voids[:,3]

        minr = radius.min()
        maxr = radius.max()
        bins = np.logspace(np.log10(minr), np.log10(maxr), nbins)
        hist, bin_edges = np.histogram(radius, bins=bins)
        r = self.bin_centres(bin_edges)

        if self.is_box:
            norm = (self.box_size**3 * np.diff(np.log10(bin_edges)))
        else:
            norm = 1

        n = hist / norm

        fig, ax = plt.subplots(1, figsize=(7,5))

        ax.plot(r, n, linestyle=linestyle, linewidth=linewidth, color=color)

        if error:
            err = np.sqrt(hist) / norm
            errplus = n + err
            errminus = n - err
            ax.fill_between(r, errplus, errminus, color='#AAAAAA')

        ax.set_xlabel(r'$R_{\rm{eff}} \hspace{0.5} /\ \rm{h^{-1}Mpc}$',
                      fontsize=17)
        ax.set_ylabel(r'$\rm{dN}\ / \rm{d} log(R_{eff})\ [h^{3}Mpc^{-3}]$',
                      fontsize=17)
        ax.set_yscale('log', nonposy='clip')
        ax.set_xlim(bins.min(), bins.max())
        ax.tick_params(which='both', width=1.0, labelsize=17)
        ax.tick_params(which='both', width=1.0, labelsize=17)

        plt.tight_layout()

        fout = self.handle + '_void_abundance.png'
        plt.savefig(fout, format='png')