
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

class CaiModel:
    '''
    Void-galaxy RSD model presented
    in Cai et al. (2016).
    '''

    def __init__(self, handle_obs, handle_mocks, mock_observation=False):

        self.handle_obs = handle_obs
        self.handle_mocks = handle_mocks

        # cosmology for Minerva
        self.om_m = 0.285
        self.s8 = 0.828
        self.cosmo = Cosmology(om_m=self.om_m, s8=self.s8)

        self.eff_z = 0.57 # effective redshift for LRGs
        self.b = 2.3 # bias for LRGs

        self.growth = self.cosmo.get_growth(self.eff_z)
        self.f = self.cosmo.get_f(self.eff_z)
        self.fs8 = self.f * self.s8 * self.growth
        self.bs8 = self.b * self.s8 * self.growth

        # build covariance matrix
        files_mocks = sorted(glob.glob(handle_mocks))
        self.s_for_xi, self.mu_for_xi, xi_smu_mocks = self.readCorrFileList(files_mocks)
        self.s_for_xi, xi2_mocks = self.getQuadrupole(self.s_for_xi, self.mu_for_xi, xi_smu_mocks)
        print(np.shape(xi2_mocks))
        cov_xi2 = self.getCovarianceMatrix(xi2_mocks)
        corr_xi2 = self.getCovarianceMatrix(xi2_mocks, norm=True)
        self.cov = cov_xi2
        self.icov = np.linalg.inv(self.cov)

        fig, ax = plt.subplots(1, figsize=(5,5))
        ax.imshow(corr_xi2, origin='lower')
        plt.savefig('/media/epaillasv/BlackIce/eboss/test/corr_xi2.png')

     
        if mock_observation:
            files_obs = sorted(glob.glob(handle_mocks))
            self.s_for_xi, self.mu_for_xi, xi_smu_obs = self.readCorrFileList(files_obs)
            xi_smu_obs = np.mean(xi_smu_obs, axis=0)
        else:
            files_obs = sorted(glob.glob(handle_obs))
            self.s_for_xi, self.mu_for_xi, xi_smu_obs = self.readCorrFile(files_obs)

        # build an interpolating bivariate spline for 
        # the correlation functions
        self.xi_smu = RectBivariateSpline(self.s_for_xi, self.mu_for_xi,
                                          xi_smu_obs, kx=3, ky=3)

        fig, ax = plt.subplots(1, figsize=(6,4))
        ax.imshow(self.xi_smu(self.s_for_xi, self.mu_for_xi))
        plt.savefig('/media/epaillasv/BlackIce/eboss/test/xi_smu.png')

        # calculate multiples for a finite set of epsilon values,
        # and build splines to use later with mcmc
        epsilon_grid = np.linspace(0.8, 1.2, 30)
        self.get_AP_splines(epsilon_grid)


        self.nwalkers = 32
        self.ndim = 1
        beta_0 = self.fs8 / self.bs8
        epsilon_0 = 1.0
        self.p0 = np.asarray([beta_0]) \
             + 1e-4*np.random.randn(self.nwalkers, self.ndim)


    def get_AP_splines(self, epsilon_grid):
        xi2_ap = []
        xi0_ap = []
        xibar_ap  = []
        alpha = 1.0
        for epsilon in epsilon_grid:
            alpha_para = alpha * epsilon ** (-2/3)
            alpha_perp = epsilon * alpha_para
            xi0, xibar, xi2 = self.theory_multipoles(alpha_perp, alpha_para,
                                                     self.s_for_xi, self.mu_for_xi)

            xi0_ap.append(xi0)
            xibar_ap.append(xibar)
            xi2_ap.append(xi2)

        xi0_ap = np.asarray(xi0_ap)
        xibar_ap = np.asarray(xibar_ap)
        xi2_ap = np.asarray(xi2_ap)

        self.xi0_APSpline = [interp1d(epsilon_grid, xi0_ap[:,i]) for i in range(len(self.s_for_xi))]
        self.xibar_APSpline = [interp1d(epsilon_grid, xibar_ap[:,i]) for i in range(len(self.s_for_xi))]
        self.xi2_APSpline = [interp1d(epsilon_grid, xi2_ap[:,i]) for i in range(len(self.s_for_xi))]

    def get_AP_multipole(self, splines, epsilon):
        xi = []
        for i in range(len(splines)):
            xi.append(splines[i](epsilon))
        return np.asarray(xi)


    def run_mcmc(self, niter=500):
        sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, self.log_probability)
        sampler.run_mcmc(self.p0, niter, progress=True);
        
        fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
        samples = sampler.get_chain()
        labels = [r"$\beta$", r"$\epsilon$"]
        for i in range(self.ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], "r", alpha=0.3)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)

        axes[-1].set_xlabel("step number")

        plt.show()

    def plot_mcmc_corner(self):
        flat_samples = self.sampler.get_chain(discard=100, thin=15, flat=True)
        fig = corner.corner(flat_samples, labels=[r"$\beta$", r"$\epsilon$"],
                    show_titles=True)

        plt.savefig('/home/epaillasv/Desktop/corner.pdf')




    def log_likelihood(self, theta):
        beta = theta
        alpha = 1.0
        epsilon = 1.0
        G = 2 * beta / (3 + beta)
        xi0 = self.get_AP_multipole(self.xi0_APSpline, epsilon)
        xibar = self.get_AP_multipole(self.xibar_APSpline, epsilon)
        xi2 = self.get_AP_multipole(self.xi2_APSpline, epsilon)

        model = G * (xi0 - xibar)
        chi2 = np.dot(np.dot((xi2 - model), self.icov), xi2 - model)
        loglike = -1/2 * chi2 - np.log((2*np.pi)**(len(self.cov)/2)) * np.sum(np.linalg.eig(self.cov)[0])
        time.sleep(5)
        print(xi2[:5])
        print(xi0[:5])
        print(np.dot((xi2 - model), self.icov))
        return loglike

    def log_prior(self, theta):
        beta = theta

        if  0.1 < beta < 1.0:
            return 0.0
        else:
            return -np.inf

    def log_probability(self, theta):
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(theta)


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
            mufunc = InterpolatedUnivariateSpline(true_mu, xi_model, k=3)
            
            # get multipoles
            monopole[i] = quad(lambda xx: mufunc(xx), 0, 1, full_output=1)[0]
            quadrupole[i] = quad(lambda xx: mufunc(xx) * 5 * (3 * xx ** 2 - 1) / 2., 0, 1, full_output=1)[0]
            
        monofunc = InterpolatedUnivariateSpline(s, monopole, k=3)
            
        # cumulative monopole
        integral = np.zeros_like(s)
        for i in range(len(integral)):
            integral[i] = quad(lambda x: monofunc(x) * x ** 2, 0, s[i], full_output=1)[0]
        monopole_bar = 3 * integral / s ** 3

        return monopole, monopole_bar, quadrupole


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
            mono[j] = quad(lambda x: mufunc(x), -1, 1, full_output=1)[0] / 2

        return s, mono

    def _getQuadrupole(self, s, mu, xi_smu):
        quadr = np.zeros(xi_smu.shape[0])
        for j in range(xi_smu.shape[0]):
            mufunc = InterpolatedUnivariateSpline(mu, xi_smu[j, :], k=3)
            quadr[j] = quad(lambda x: mufunc(x) * 5 * (3. * x ** 2 - 1) / 2., -1, 1, full_output=1)[0] / 2

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
            



