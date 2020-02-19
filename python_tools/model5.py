
import numpy as np
import sys
import os
import glob
import subprocess
from astropy.io import fits
from python_tools.cosmology import Cosmology
from python_tools.galaxycat import GalaxyCatalogue
from scipy.spatial import Delaunay
from scipy.integrate import quad
from scipy.interpolate import RectBivariateSpline, InterpolatedUnivariateSpline, interp1d
import emcee
import corner
from scipy.io import FortranFile
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import time




class Model5:
    '''
    Void-galaxy RSD model presented
    in Cai et al. (2016).
    '''

    def __init__(self, delta_r_file, xi_r_file, sv_file, xi_smu_file,
                 covmat_file, xi_smu_mocks=''):

        self.delta_r_file = delta_r_file
        self.xi_r_file = xi_r_file
        self.sv_file = sv_file
        self.xi_smu_file = xi_smu_file
        self.covmat_file = covmat_file
        self.xi_smu_mocks = xi_smu_mocks


        print("Setting up Void RSD model #5 .")

        # cosmology for Minerva
        self.om_m = 0.285
        self.s8 = 0.828
        self.cosmo = Cosmology(om_m=self.om_m, s8=self.s8)

        self.eff_z = 0.57
        self.b = 2.01 

        self.growth = self.cosmo.get_growth(self.eff_z)
        self.f = self.cosmo.get_f(self.eff_z)
        self.s8norm = self.s8 * self.growth 

        eofz = np.sqrt((self.om_m * (1 + self.eff_z) ** 3 + 1 - self.om_m))
        self.iaH = (1 + self.eff_z) / (100. * eofz) 

        # build covariance matrix
        if os.path.isfile(self.covmat_file):
            print('Reading covariance matrix: ' + self.covmat_file)
            self.cov = np.load(self.covmat_file)
        else:
            print('Computing covariance matrix...')
            self.cov = self.MultipoleCovariance()
            np.save(self.covmat_file, self.cov)

        self.icov = np.linalg.inv(self.cov)

        # read real-space monopole
        data = np.genfromtxt(self.xi_r_file)
        self.r_for_xi = data[:,0]
        xi_r = data[:,-1]
        self.xi_r = InterpolatedUnivariateSpline(self.r_for_xi, xi_r, k=1, ext=3)

        # read void-matter correlation function
        data = np.genfromtxt(self.delta_r_file)
        self.r_for_delta = data[:,0]
        delta_r = data[:,-1]
        self.delta_r = InterpolatedUnivariateSpline(self.r_for_delta, delta_r, k=1, ext=3)

        integral = np.zeros_like(self.r_for_delta)
        for i in range(len(integral)):
            integral[i] = quad(lambda x: self.delta_r(x) * x ** 2, 0, self.r_for_delta[i], full_output=1)[0]
        Delta_r = 3 * integral / self.r_for_delta ** 3
        self.Delta_r = InterpolatedUnivariateSpline(self.r_for_delta, Delta_r, k=1, ext=3)

        # read los velocity dispersion profile
        data = np.genfromtxt(self.sv_file)
        self.r_for_sv = data[:,0]
        sv = data[:,-1] / data[-1, -1]
        self.sv = InterpolatedUnivariateSpline(self.r_for_sv, sv, k=1, ext=3)

        # read redshift-space correlation function
        s, mu, xi_smu_obs = self.readCorrFile(self.xi_smu_file)
        s, self.xi0_s = self._getMonopole(s, mu, xi_smu_obs)
        s, self.xi2_s = self._getQuadrupole(s, mu, xi_smu_obs)

        self.datavec = np.concatenate((self.xi0_s, self.xi2_s))

        self.s_for_xi = s
        self.mu_for_xi = mu

    
    def log_likelihood(self, theta):
        fs8, bs8, sigma_v, epsilon = theta
        alpha = 1.0
        alpha_para = alpha * epsilon ** (-2/3)
        alpha_perp = epsilon * alpha_para

        xi0, xibar, xi2 = self.theory_multipoles(fs8, bs8, sigma_v,
                                                 alpha_perp, alpha_para,
                                                 self.s_for_xi, self.mu_for_xi)

        datavec = np.concatenate((xi0, xi2))

        chi2 = np.dot(np.dot((self.datavec - datavec), self.icov), self.datavec - datavec)
        loglike = -1/2 * chi2 - np.log((2*np.pi)**(len(self.cov)/2)) * np.sum(np.linalg.eig(self.cov)[0])
        return loglike

    def log_prior(self, theta):
        fs8, bs8, sigma_v, epsilon = theta
        beta = fs8 / bs8

        if 0.1 < fs8 < 0.8 and 0.45 < bs8 < 1.9 and 250 < sigma_v < 500 \
        and 0.8 < epsilon < 1.2:
            return 0.0
        
        return -np.inf


    def theory_multipoles(self, fs8, bs8, sigma_v, alpha_perp, alpha_para, s, mu):

        monopole = np.zeros(len(s))
        quadrupole = np.zeros(len(s))
        true_mu = np.zeros(len(mu))
        xi_model = np.zeros(len(mu))
        scaled_fs8 = fs8 / self.s8norm

        # rescale input monopole functions to account for alpha values
        mus = np.linspace(0, 1., 101)
        r = self.r_for_delta
        rescaled_r = np.zeros_like(r)
        for i in range(len(r)):
            rescaled_r[i] = np.trapz((r[i] * alpha_para) * np.sqrt(1. + (1. - mus ** 2) *
                            (alpha_perp ** 2 / alpha_para ** 2 - 1)), mus)

        x = rescaled_r
        y1 = self.xi_r(r)
        y2 = self.delta_r(r)
        y3 = self.Delta_r(r)
        y4 = self.sv(r)

        # build rescaled interpolating functions using the relabelled separation vectors
        rescaled_xi_r = InterpolatedUnivariateSpline(x, y1, k=3)
        rescaled_delta_r = InterpolatedUnivariateSpline(x, y2, k=3, ext=3)
        rescaled_Delta_r = InterpolatedUnivariateSpline(x, y3, k=3, ext=3)
        rescaled_sv = InterpolatedUnivariateSpline(x, y4, k=3, ext=3)
        sigma_v = alpha_para * sigma_v
 

        for i in range(len(s)):
            for j in range(len(mu)):
                true_sperp = s[i] * np.sqrt(1 - mu[j] ** 2) * alpha_perp
                true_spar = s[i] * mu[j] * alpha_para
                true_s = np.sqrt(true_spar ** 2. + true_sperp ** 2.)
                true_mu[j] = true_spar / true_s

                rpar = true_spar + true_s * scaled_fs8 * rescaled_Delta_r(true_s) * true_mu[j] / 3.
                sy_central = sigma_v * rescaled_sv(np.sqrt(true_sperp**2 + rpar**2)) * self.iaH
                y = np.linspace(-3 * sy_central, 3 * sy_central, 100)

                rpar = true_spar + true_s * scaled_fs8 * rescaled_Delta_r(true_s) * true_mu[j] / 3. - y
                rr = np.sqrt(true_sperp ** 2 + rpar ** 2)
                sy = sigma_v * rescaled_sv(rr) * self.iaH

                integrand = (1 + rescaled_xi_r(rr)) * \
                            (1 + (scaled_fs8 * rescaled_Delta_r(rr) / 3. - y * true_mu[j] / rr) * (1 - true_mu[j]**2) +
                             scaled_fs8 * (rescaled_delta_r(rr) - 2 * rescaled_Delta_r(rr) / 3.) * true_mu[j]**2)
                integrand = integrand * np.exp(-(y**2) / (2 * sy**2)) / (np.sqrt(2 * np.pi) * sy)
                xi_model[j] = np.trapz(integrand, y) - 1


            # build interpolating function for xi_smu at true_mu
            mufunc = InterpolatedUnivariateSpline(true_mu, xi_model, k=3)
            
            # get multipoles
            monopole[i] = quad(lambda xx: mufunc(xx) / 2, -1, 1, full_output=1)[0]
            quadrupole[i] = quad(lambda xx: mufunc(xx) * 5 / 2* (3 * xx ** 2 - 1) / 2., -1, 1, full_output=1)[0]
            
        monofunc = InterpolatedUnivariateSpline(s, monopole, k=3)
            
        # cumulative monopole
        integral = np.zeros_like(s)
        for i in range(len(integral)):
            integral[i] = quad(lambda x: monofunc(x) * x ** 2, 0, s[i], full_output=1)[0]
        monopole_bar = 3 * integral / s ** 3

        return monopole, monopole_bar, quadrupole

    
    def MultipoleCovariance(self):
        if self.xi_smu_mocks == '':
            sys.exit('Could not find mock files to build covariance. Aborting...')

        files_mocks = sorted(glob.glob(self.xi_smu_mocks))
        mock_datavec = []
        for fname in files_mocks:
            s, mu, xi_smu_mock = self.readCorrFile(fname)
            s, xi0 = self._getMonopole(s, mu, xi_smu_mock)
            s, xi2 = self._getQuadrupole(s, mu, xi_smu_mock)

            datavec = np.concatenate((xi0, xi2))

            mock_datavec.append(datavec)

        mock_datavec = np.asarray(mock_datavec)
        cov = self.CovarianceMatrix(mock_datavec)
        
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


            



