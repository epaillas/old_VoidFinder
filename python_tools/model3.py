
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

class Model3:
    '''
    Void-galaxy RSD model presented
    in Cai et al. (2016).
    '''

    def __init__(self, xi_smu_file, covmat_file, xi_smu_mocks=''):

        self.xi_smu_file = xi_smu_file
        self.xi_smu_mocks = xi_smu_mocks
        self.covmat_file = covmat_file
        self.nmocks = 120

        print("Setting up void RSD model #2.")

        # cosmology for Minerva
        self.om_m = 0.285
        self.s8 = 0.828
        self.cosmo = Cosmology(om_m=self.om_m, s8=self.s8)

        self.eff_z = 0.57 # effective redshift for LRGs
        self.b = 2.05 # bias for LRGs

        self.growth = self.cosmo.get_growth(self.eff_z)
        self.f = self.cosmo.get_f(self.eff_z)
        self.fs8 = self.f * self.s8 * self.growth
        self.bs8 = self.b * self.s8 * self.growth
        self.beta = self.fs8 / self.bs8
        self.G = (2 * self.beta) / (3 + self.beta)

        eofz = np.sqrt((self.om_m * (1 + self.eff_z) ** 3 + 1 - self.om_m))
        self.iaH = (1 + self.eff_z) / (100. * eofz)

        print('1/aH = {}'.format(self.iaH))
        print('f = {}'.format(self.f))
        print('beta = {}'.format(self.beta))
        print('G = {}'.format(self.G))
        print('fs8 = {}'.format(self.fs8))
        print('bs8 = {}'.format(self.bs8))
        print('growth = {}'.format(self.growth))

        # build covariance matrix
        if os.path.isfile(self.covmat_file):
            print('Reading covariance matrix: ' + self.covmat_file)
            self.cov = np.load(self.covmat_file)
        else:
            print('Computing covariance matrix...')
            self.cov = self.MultipoleCovariance()
            np.save(self.covmat_file, self.cov)

        self.icov = np.linalg.inv(self.cov)

        # build an interpolating spline for xi_smu
        self.s_for_xi, self.mu_for_xi, xi_smu_obs = self.readCorrFile(self.xi_smu_file)
        self.xi_smu = RectBivariateSpline(self.s_for_xi, self.mu_for_xi,
                                          xi_smu_obs, kx=3, ky=3)



    def log_likelihood(self, theta):
        beta, epsilon = theta
        alpha = 1.0
        alpha_para = alpha * epsilon ** (-2/3)
        alpha_perp = epsilon * alpha_para
        G = 2 * beta / (3 + beta)

        xi0, xibar, xi2 = self.theory_multipoles(alpha_perp, alpha_para,
                                                 self.s_for_xi, self.mu_for_xi)

        model = G * (xi0 - xibar)
        chi2 = np.dot(np.dot((xi2 - model), self.icov), xi2 - model)
        loglike = -self.nmocks/2 * np.log(1 + chi2/(self.nmocks-1))
        return loglike

    def log_prior(self, theta):
        beta, epsilon = theta

        if  0.1 < beta < 1.0 and 0.8 < epsilon < 1.2:
            return 0.0
        else:
            return -np.inf


    def theory_multipoles(self, alpha_perp, alpha_para, s, mu):

        monopole = np.zeros(len(s))
        quadrupole = np.zeros(len(s))
        true_mu = np.zeros(len(mu))
        xi_model = np.zeros(len(mu))

        for i in range(len(s)):
            for j in range(len(mu)):
                true_sperp = s[i] * np.sqrt(1 - mu[j] ** 2) * alpha_perp
                true_spar = s[i] * mu[j] * alpha_para
                true_s = np.sqrt(true_spar ** 2. + true_sperp ** 2.)
                true_mu[j] = true_spar / true_s

                xi_model[j] = self.xi_smu(true_s, true_mu[j])


            # build interpolating function for xi_smu at true_mu
            mufunc = InterpolatedUnivariateSpline(true_mu[np.argsort(true_mu)],
                                                  xi_model[np.argsort(true_mu)], k=3)

            # get multipoles
            xaxis = np.linspace(-1, 1, 1000)

            yaxis = mufunc(xaxis) / 2
            monopole[i] = simps(yaxis, xaxis)

            yaxis = mufunc(xaxis) * 5 / 2 * (3 * xaxis**2 - 1) / 2
            quadrupole[i] = simps(yaxis, xaxis)
            
            
        monofunc = InterpolatedUnivariateSpline(s, monopole, k=3)
            
        # cumulative monopole
        integral = np.zeros_like(s)
        for i in range(len(integral)):
            xaxis = np.linspace(0, s[i], 1000)
            yaxis = monofunc(xaxis) * xaxis**2
            integral[i] = simps(yaxis, xaxis)
        monopole_bar = 3 * integral / s ** 3

        return monopole, monopole_bar, quadrupole

    def MultipoleCovariance(self):
        files_mocks = sorted(glob.glob(self.xi_smu_mocks))
        mock_xi0 = []
        mock_xi2 = []
        for fname in files_mocks:
            s, mu, xi_smu_mock = self.readCorrFile(fname)
            s, xi0 = self._getMonopole(s, mu, xi_smu_mock)
            s, xi2 = self._getQuadrupole(s, mu, xi_smu_mock)

            mock_xi0.append(xi0)
            mock_xi2.append(xi2)

        mock_xi0 = np.asarray(mock_xi0)
        mock_xi2 = np.asarray(mock_xi2)

        cov_xi0 = self.CovarianceMatrix(mock_xi0)
        cov_xi2 = self.CovarianceMatrix(mock_xi2)
        cov_xi02 = self.CrossCovarianceMatrix(mock_xi0, mock_xi2)
        cov_xi20 = self.CrossCovarianceMatrix(mock_xi2, mock_xi0)

        cov = cov_xi2 + self.G**2 * cov_xi0 -\
              self.G * cov_xi02 - self.G * cov_xi20

        return cov


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

    def CovarianceMatrix(self, data, norm=False):
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

    def CrossCovarianceMatrix(self, data1, data2, norm=False):
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



            



