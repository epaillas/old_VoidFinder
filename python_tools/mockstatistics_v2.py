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
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.integrate import quad

class MockStatistics:

    def __init__(self, handle):

        self.handle = handle

        self.bin_sampling = 2
        self.cutoff_lo = 1
        self.cutoff_hi = None

        self.files = glob.glob(handle)
        self.data = self.getData()
        self.bins, self.monopole = self.getMonopole()
        self.bins, self.quadrupole = self.getQuadrupole()
        self.bins, self.average_monopole = self.getAverageMonopole()
        self.bins, self.differential_monopole = self.getDiffMonopole()

        
    def getData(self):
        data = []
        for fname in self.files[:5]:
            f = np.genfromtxt(fname)
            f[np.isnan(f)] = -1 # Set invalid values to -1
            f[f == np.inf] = -1 # Set invalid values to -1     
            data.append(f)
        data = np.asarray(data)
        print(np.shape(data))
        return data

    def gluedBins(self, bins, factor):
        glued_bins = []
        for i in range(0, len(bins), factor):
            glued_bins.append(np.mean(bins[i:i+factor]))
        return np.asarray(glued_bins)

    def qFactor(self, mu, qper, qpara):
        q = np.sqrt(qpara**2 * mu**2 + qper**2*(1 - mu**2))
        return q

    def _getMonopole(self, data, epsilon=1):
        '''
        Computes the monopole from an input
        3D correlation function. Distorts
        input s and mu vectors according to the
        value of the epsilon AP parameter.
        '''

        s = np.unique(data[:,0])
        mu = np.unique(data[:,3])
            
        qpara = 1 * epsilon ** (-2/3)
        qperp = epsilon * qpara

        true_sperp = s * np.sqrt(1 - mu**2) * qperp
        true_spara = s * mu * qpara
        true_s = np.sqrt(true_spara ** 2 + true_sperp ** 2)
        true_mu = true_spara / true_s

        xi_smu = np.zeros([len(s), len(mu)])
        counter = 0
        for i in range(len(s)):
            for j in range(len(mu)):
                xi_smu[i, j] = data[counter, -1]
                counter += 1

        mono = np.zeros(xi_smu.shape[0])
        for j in range(xi_smu.shape[0]):
            mufunc = InterpolatedUnivariateSpline(true_mu, xi_smu[j, :], k=3)
            mono[j] = quad(lambda x: mufunc(x), 0, 1, full_output=1)[0]


        return true_s, mono

    def _getQuadrupole(self, data, epsilon=1):
        s = np.unique(data[:,0])
        mu = np.unique(data[:,3])

        qpara = 1 * epsilon ** (-2/3)
        qperp = epsilon * qpara

        true_sperp = s * np.sqrt(1 - mu**2) * qperp
        true_spara = s * mu * qpara
        true_s = np.sqrt(true_spara ** 2 + true_sperp ** 2)
        true_mu = true_spara / true_s

        xi_smu = np.zeros([len(s), len(mu)])
        counter = 0
        for i in range(len(s)):
            for j in range(len(mu)):
                xi_smu[i, j] = data[counter, -1]
                counter += 1

        quadr = np.zeros(xi_smu.shape[0])
        for j in range(xi_smu.shape[0]):
            mufunc = InterpolatedUnivariateSpline(true_mu, xi_smu[j, :], k=3)
            quadr[j] = quad(lambda x: mufunc(x) * 5 * (3. * x ** 2 - 1) / 2., 0, 1, full_output=1)[0]

        return true_s, quadr

    def _getAverageMonopole(self, data, epsilon=1):

        s = np.unique(data[:,0])
        mu = np.unique(data[:,3])

        qpara = 1 * epsilon ** (-2/3)
        qperp = epsilon * qpara

        true_sperp = s * np.sqrt(1 - mu**2) * qperp
        true_spara = s * mu * qpara
        true_s = np.sqrt(true_spara ** 2 + true_sperp ** 2)
        true_mu = true_spara / true_s

        xi_smu = np.zeros([len(s), len(mu)])
        counter = 0
        for i in range(len(s)):
            for j in range(len(mu)):
                xi_smu[i, j] = data[counter, -1]
                counter += 1

        mono = np.zeros(xi_smu.shape[0])
        for j in range(xi_smu.shape[0]):
            mufunc = InterpolatedUnivariateSpline(true_mu, xi_smu[j, :], k=3)
            mono[j] = quad(lambda x: mufunc(x), 0, 1, full_output=1)[0]

        delta_r = InterpolatedUnivariateSpline(true_s, mono, k=3)

        Delta_r = np.zeros(xi_smu.shape[0])
        for i in range(len(true_s)):
            integral = quad(lambda x: delta_r(x) * x ** 2, 0, true_s[i], full_output=1)[0]
            Delta_r[i] = (3 / true_s[i]**3) * integral

        return true_s, Delta_r

    def getMonopole(self, epsilon=1):
        monopole = []
        for data in self.data:
            r, xi0 = self._getMonopole(data, epsilon=epsilon)
            r = self.gluedBins(r, self.bin_sampling)[self.cutoff_lo:self.cutoff_hi]
            xi0 = self.gluedBins(xi0, self.bin_sampling)[self.cutoff_lo:self.cutoff_hi]
            monopole.append(xi0)
        monopole = np.asarray(monopole)
        return r, monopole

    def getQuadrupole(self, epsilon=1):
        quadrupole = []
        for data in self.data:
            r, xi2 = self._getQuadrupole(data, epsilon=epsilon)
            r = self.gluedBins(r, self.bin_sampling)[self.cutoff_lo:self.cutoff_hi]
            xi2 = self.gluedBins(xi2, self.bin_sampling)[self.cutoff_lo:self.cutoff_hi]
            quadrupole.append(xi2)
        quadrupole = np.asarray(quadrupole)
        return r, quadrupole

    def getAverageMonopole(self, epsilon=1):
        avmono = []
        for data in self.data:
            r, xiav = self._getAverageMonopole(data, epsilon=epsilon)
            r = self.gluedBins(r, self.bin_sampling)[self.cutoff_lo:self.cutoff_hi]
            xiav = self.gluedBins(xiav, self.bin_sampling)[self.cutoff_lo:self.cutoff_hi]
            avmono.append(xiav)
        avmono = np.asarray(avmono)
        return r, avmono

    def getDiffMonopole(self, epsilon=1):
        diffmono = []
        for data in self.data:
            r, xi0 = self._getMonopole(data, epsilon=epsilon)
            r, xiav = self._getAverageMonopole(data, epsilon=epsilon)
            dxi0 = xi0 - xiav
            r = self.gluedBins(r, self.bin_sampling)[self.cutoff_lo:self.cutoff_hi]
            dxi0 = self.gluedBins(dxi0, self.bin_sampling)[self.cutoff_lo:self.cutoff_hi]
            diffmono.append(dxi0)
        diffmono = np.asarray(diffmono)
        return r, diffmono

    def getCovarianceMatrix(self, data, norm=False):
        """
        Assumes rows are observations,
        columns are variables
        """
        nobs, nbins = np.shape(data)
        mean = np.mean(data, axis=0)
        cov = np.zeros([nbins, nbins])

        for k in range(nobs):
            for i in range(nbins):
                for j in range(nbins):
                    cov[i, j] += (data[k, i] - mean[i])*(data[k, j] - mean[j])

        cov /= nobs - 1
        
        if norm:
            corr = np.zeros_like(cov)
            for i in range(nbins):
                for j in range(nbins):
                    corr[i, j] = cov[i, j] / np.sqrt(cov[i, i] * cov[j, j])
            return corr
        else:
            return cov

    def getCrossCovarianceMatrix(self, data1, data2, norm=False):
        """
        Assumes rows are observations,
        columns are variables
        """
        nobs, nbins = np.shape(data1)
        mean1 = np.mean(data1, axis=0)
        mean2 = np.mean(data2, axis=0)
        cov = np.zeros([nbins, nbins])

        for k in range(nobs):
            for i in range(nbins):
                for j in range(nbins):
                    cov[i, j] += (data1[k, i] - mean1[i])*(data2[k, j] - mean2[j])

        cov /= nobs - 1
        
        if norm:
            corr = np.zeros_like(cov)
            for i in range(nbins):
                for j in range(nbins):
                    corr[i, j] = cov[i, j] / np.sqrt(cov[i, i] * cov[j, j])
            return corr
        else:
            return cov
            



