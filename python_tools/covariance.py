import numpy as np
import sys
import os
import glob
import subprocess
from multiprocessing import Pool, cpu_count
from astropy.io import fits
from python_tools.cosmology import Cosmology
from python_tools.galaxycat import GalaxyCatalogue
from scipy.spatial import Delaunay
from scipy.integrate import quad, simps
from scipy.interpolate import RectBivariateSpline, InterpolatedUnivariateSpline, interp1d
import emcee
import corner
from scipy.io import FortranFile
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import time

def readCorrFileList(self, fnames):
        xi_smu_list = []
        for fname in fnames:
            data = np.genfromtxt(fname)
            data[np.isnan(data)] = -1 # Set invalid values to -1
            data[data == np.inf] = -1 # Set invalid values to -1 
            s = np.unique(data[:,0])
            mu = np.unique(data[:,3])

            xi_smu = np.zeros([len(s), len(mu)])
            counter = 0
            for i in range(len(s)):
                for j in range(len(mu)):
                    xi_smu[i, j] = data[counter, -1]
                    counter += 1

            xi_smu_list.append(xi_smu)

        xi_smu_list = np.asarray(xi_smu_list)
        return s, mu, xi_smu_list

    def readCorrFile(self, fname):
        data = np.genfromtxt(fname)
        s = np.unique(data[:,0])
        mu = np.unique(data[:,3])

        xi_smu = np.zeros([len(s), len(mu)])
        counter = 0
        for i in range(len(s)):
            for j in range(len(mu)):
                xi_smu[i, j] = data[counter, -1]
                counter += 1

        return s, mu, xi_smu

    def getMonopole(self, s, mu, xi_smu_list):
        monopole = []
        for xi_smu in xi_smu_list:
            s, xi0 = self._getMonopole(s, mu, xi_smu)
            monopole.append(xi0)
        monopole = np.asarray(monopole)
        return s, monopole

    def getQuadrupole(self, s, mu, xi_smu_list):
        quadrupole = []
        for xi_smu in xi_smu_list:
            r, xi2 = self._getQuadrupole(s, mu, xi_smu)
            quadrupole.append(xi2)
        quadrupole = np.asarray(quadrupole)
        return r, quadrupole

    def _getMonopole(self, s, mu, xi_smu):
        mono = np.zeros(xi_smu.shape[0])
        for j in range(xi_smu.shape[0]):
            mufunc = InterpolatedUnivariateSpline(mu, xi_smu[j, :], k=3)
            mono[j] = quad(lambda x: mufunc(x) / 2, -1, 1, full_output=1)[0]

        return s, mono

    def _getQuadrupole(self, s, mu, xi_smu):
        quadr = np.zeros(xi_smu.shape[0])
        for j in range(xi_smu.shape[0]):
            mufunc = InterpolatedUnivariateSpline(mu, xi_smu[j, :], k=3)
            quadr[j] = quad(lambda x: mufunc(x) * 5 / 2 * (3. * x ** 2 - 1) / 2., -1, 1, full_output=1)[0]

        return s, quadr

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

def MultipoleCovariance_Model1(handle_mocks, handle_cov):
        files_mocks = sorted(glob.glob(handle_mocks))
        mock_xi0 = []
        mock_dxi0 = []
        mock_xi2 = []
        for fname in files_mocks:
            s, mu, xi_smu_mock = readCorrFile(fname)
            s, xi0 = _getMonopole(s, mu, xi_smu_mock)
            s, xi2 = _getQuadrupole(s, mu, xi_smu_mock)

            monofunc = InterpolatedUnivariateSpline(s, xi0, k=3)
            integral = np.zeros_like(s)
            for i in range(len(integral)):
                integral[i] = quad(lambda x: monofunc(x) * x ** 2, 0, s[i], full_output=1)[0]
            xibar = 3 * integral / s ** 3

            dxi0 = xi0 - xibar

            mock_xi0.append(xi0)
            mock_dxi0.append(dxi0)
            mock_xi2.append(xi2)

        mock_xi0 = np.asarray(mock_xi0)
        mock_dxi0 = np.asarray(mock_dxi0)
        mock_xi2 = np.asarray(mock_xi2)

        cov_xi0 = getCovarianceMatrix(mock_xi0)
        cov_dxi0 = getCovarianceMatrix(mock_dxi0)
        cov_xi2 = getCovarianceMatrix(mock_xi2)
        cov_xi02 = getCrossCovarianceMatrix(mock_dxi0, mock_xi2)
        cov_xi20 = getCrossCovarianceMatrix(mock_xi2, mock_dxi0)

        cout = [cov_xi0, cov_dxi0, cov_xi2, cov_xi02, cov_xi20]

        return cov_xi0, cov_dxi0, cov_xi2, cov_xi02, cov_xi20

def MultipoleCovariance_Model5(handle_mocks, handle_cov):
        files_mocks = sorted(glob.glob(handle_mocks))
        mock_datavec = []
        for fname in files_mocks:
            s, mu, xi_smu_mock = readCorrFile(fname)
            s, xi0 = _getMonopole(s, mu, xi_smu_mock)
            s, xi2 = _getQuadrupole(s, mu, xi_smu_mock)

            datavec = np.concatenate(xi0, xi2)

            mock_datavec.append(datavec)

        mock_datavec = np.asarray(mock_datavec)

        cov = getCovarianceMatrix(mock_datavec)

        return cov

handle_mocks = '/Volumes/BlackIce/eboss/minerva_cats/void_stats/\
real_real/\
Galaxies_HOD_*_z0.57_Real_Real.SVF_recen_ovl0.5_MedianVoids.VG_CCF_rmu'

cov_model1 = MultipoleCovariance_Model1(handle_mocks=handle_mocks)

np.save(fout, cov_model1)