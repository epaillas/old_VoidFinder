
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

    def __init__(self, handle_obs, handle_mocks):

        self.handle_obs = handle_obs
        self.handle_mocks = handle_mocks

        print("Setting up Cai's void RSD model.")
        print('handle_obs: ' + handle_obs)
        print('handle_mocks: ' + handle_mocks)

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
        self.beta = self.fs8 / self.bs8
        self.G = (2 * self.beta) / (3 + self.beta)

        # build covariance matrix
        self.handle_cov = self.handle_obs + '_covmat.npy'
        if os.path.isfile(self.handle_cov):
            cov_list = np.load(self.handle_cov)
        
        else:
            cov_list = self.getMultipoleCovariance()

        cov_dxi0, cov_xi2, cov_xi02, cov_xi20 = cov_list

        self.cov = cov_xi2 + self.G**2 * cov_dxi0 - self.G*cov_xi02 - self.G*cov_xi20
        self.icov = np.linalg.inv(self.cov)


        # build an interpolating bivariate spline for 
        # the correlation functions
        self.s_for_xi, self.mu_for_xi, xi_smu_obs = self.readCorrFile(handle_obs)
        self.xi_smu = RectBivariateSpline(self.s_for_xi, self.mu_for_xi,
                                          xi_smu_obs, kx=3, ky=3)


        # calculate multiples for a finite set of epsilon values,
        # and build splines to use later with mcmc
        #epsilon_grid = np.linspace(0.8, 1.2, 30)
        #self.get_AP_splines(epsilon_grid)


        #self.nwalkers = 64
        #self.ndim = 2
        #beta_0 = self.fs8 / self.bs8
        #epsilon_0 = 1.0
        #self.p0 = np.asarray([beta_0, epsilon_0]) \
        #     + 1e-4*np.random.randn(self.nwalkers, self.ndim)

    def getMultipoleCovariance(self):
        files_mocks = sorted(glob.glob(handle_mocks))
        mock_dxi0 = []
        mock_xi2 = []
        for fname in files_mocks:
            s, mu, xi_smu_mock = self.readCorrFile(fname)
            s, xi0 = self._getMonopole(s, mu, xi_smu_mock)
            s, xi2 = self._getQuadrupole(s, mu, xi_smu_mock)

            monofunc = InterpolatedUnivariateSpline(s, xi0, k=3)
            integral = np.zeros_like(s)
            for i in range(len(integral)):
                integral[i] = quad(lambda x: monofunc(x) * x ** 2, 0, s[i], full_output=1)[0]
            xibar = 3 * integral / s ** 3

            dxi0 = xi0 - xibar

            mock_dxi0.append(dxi0)
            mock_xi2.append(xi2)

        mock_dxi0 = np.asarray(mock_dxi0)
        mock_xi2 = np.asarray(mock_xi2)

        cov_dxi0 = self.getCovarianceMatrix(mock_dxi0)
        cov_xi2 = self.getCovarianceMatrix(mock_xi2)
        cov_xi02 = self.getCrossCovarianceMatrix(mock_dxi0, mock_xi2)
        cov_xi20 = self.getCrossCovarianceMatrix(mock_xi2, mock_dxi0)

        cout = [cov_dxi0, cov_xi2, cov_xi02, cov_xi20]
        np.save(self.handle_cov, cout)

        return cov_dxi0, cov_xi2, cov_xi02, cov_xi20

    def run_mcmc(self, niter=500, backend_name=''):
        if backend_name == '':
            backend_name = self.handle_obs + '_emceeChain.h5'

        print('Running emcee with the following parameters:')
        print('nwalkers: ' + str(self.nwalkers))
        print('ndim: ' + str(self.ndim))
        print('niter: ' + str(niter))
        print('backend: ' + backend_name)


        backend = emcee.backends.HDFBackend(backend_name)
        backend.reset(self.nwalkers, self.ndim)

        self.sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim,
                                             self.log_probability,
                                             backend=backend)
        self.sampler.run_mcmc(self.p0, niter, progress=True)

    def plot_mcmc_stats(self):

        print('Saving emcee statistics...')

        fig, axes = plt.subplots(self.ndim, figsize=(10, 3.5*self.ndim),
                                 sharex=True)
        samples = self.sampler.get_chain()
        labels = [r"$\beta$", r"$\epsilon$"]
        for i in range(self.ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], "r", alpha=0.3)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)

        axes[-1].set_xlabel("Step number")

        fout = self.handle_obs + '_emceeFlatChains.png'
        print('Saving chains: ' + fout)
        plt.savefig(fout)

        flat_samples = self.sampler.get_chain(discard=100, thin=15, flat=True)
        fig = corner.corner(flat_samples, labels=[r"$\beta$", r"$\epsilon$"],
                    show_titles=True, quantiles=[0.16, 0.84],
                    truths=[self.fs8/self.bs8, 1.0], truth_color='r')
        fout = self.handle_obs + '_emceeCorner.png'
        print('Saving corner: ' + fout)
        plt.savefig(fout)




    def log_likelihood(self, theta):
        beta, epsilon = theta
        alpha = 1.0
        G = 2 * beta / (3 + beta)
        xi0 = self.get_AP_multipole(self.xi0_APSpline, epsilon)
        xibar = self.get_AP_multipole(self.xibar_APSpline, epsilon)
        xi2 = self.get_AP_multipole(self.xi2_APSpline, epsilon)

        model = G * (xi0 - xibar)
        chi2 = np.dot(np.dot((xi2 - model), self.icov), xi2 - model)
        loglike = -1/2 * chi2 - np.log((2*np.pi)**(len(self.cov)/2)) * np.sum(np.linalg.eig(self.cov)[0])
        return loglike

    def log_prior(self, theta):
        beta, epsilon = theta

        if  0.1 < beta < 1.0 and 0.8 < epsilon < 1.2:
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
            monopole[i] = quad(lambda xx: mufunc(xx) / 2, -1, 1, full_output=1)[0]
            quadrupole[i] = quad(lambda xx: mufunc(xx) * 5 / 2* (3 * xx ** 2 - 1) / 2., -1, 1, full_output=1)[0]
            
        monofunc = InterpolatedUnivariateSpline(s, monopole, k=3)
            
        # cumulative monopole
        integral = np.zeros_like(s)
        for i in range(len(integral)):
            integral[i] = quad(lambda x: monofunc(x) * x ** 2, 0, s[i], full_output=1)[0]
        monopole_bar = 3 * integral / s ** 3

        return monopole, monopole_bar, quadrupole

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






class NadathurModel:
    '''
    Void-galaxy RSD model presented
    in Cai et al. (2016).
    '''

    def __init__(self, delta_r_file, xi_r_file, sv_file, xi_smu_file):

        self.delta_r_file = delta_r_file
        self.xi_r_file = xi_r_file
        self.sv_file = sv_file
        self.xi_smu_file = xi_smu_file


        print("Setting up Nadathur's void RSD model.")

        # cosmology for Minerva
        self.om_m = 0.285
        self.s8 = 0.828
        self.cosmo = Cosmology(om_m=self.om_m, s8=self.s8)

        self.eff_z = 0.57 # effective redshift for LRGs
        self.b = 2.3 # bias for LRGs

        self.growth = self.cosmo.get_growth(self.eff_z)
        self.f = self.cosmo.get_f(self.eff_z)
        self.s8norm = self.s8 * self.growth 

        # build covariance matrix
        # self.handle_cov = self.handle_real_redshift_mocks + '_covmat.npy'
        # if os.path.isfile(self.handle_cov):
        #     cov_list = np.load(self.handle_cov)
        
        # else:
        #     cov_list = self.getMultipoleCovariance()

        # cov_dxi0, cov_xi2, cov_xi02, cov_xi20 = cov_list

        # self.cov = cov_xi2 + self.G**2 * cov_dxi0 - self.G*cov_xi02 - self.G*cov_xi20
        # self.icov = np.linalg.inv(self.cov)


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
        sv = data[:,-1] / data[-1, 1]
        self.sv = InterpolatedUnivariateSpline(self.r_for_sv, sv, k=1, ext=3)

        # read redshift-space correlation function
        s, mu xi_smu_obs = self.readCorrFile(self.xi_smu_file)
        s, self.xi0 = self._getMonopole(s, mu, xi_smu_obs)
        s, self.xi2 = self._getQuadrupole(s, mu, xi_smu_obs)

        monofunc = InterpolatedUnivariateSpline(s, xi0, k=3)
        integral = np.zeros_like(s)
        for i in range(len(integral)):
            integral[i] = quad(lambda x: monofunc(x) * x ** 2, 0, s[i], full_output=1)[0]
        self.xibar = 3 * integral / s ** 3

        # build an interpolating bivariate spline for 
        # the correlation functions
        # self.xi_smu = RectBivariateSpline(self.s_for_xi, self.mu_for_xi,
        #                                   xi_smu_obs, kx=3, ky=3)



        #self.nwalkers = 64
        #self.ndim = 2
        #beta_0 = self.fs8 / self.bs8
        #epsilon_0 = 1.0
        #self.p0 = np.asarray([beta_0, epsilon_0]) \
        #     + 1e-4*np.random.randn(self.nwalkers, self.ndim)

    

    def run_mcmc(self, niter=500, backend_name=''):
        if backend_name == '':
            backend_name = self.handle_obs + '_emceeChain.h5'

        print('Running emcee with the following parameters:')
        print('nwalkers: ' + str(self.nwalkers))
        print('ndim: ' + str(self.ndim))
        print('niter: ' + str(niter))
        print('backend: ' + backend_name)


        backend = emcee.backends.HDFBackend(backend_name)
        backend.reset(self.nwalkers, self.ndim)

        self.sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim,
                                             self.log_probability,
                                             backend=backend)
        self.sampler.run_mcmc(self.p0, niter, progress=True)

    def plot_mcmc_stats(self):

        print('Saving emcee statistics...')

        fig, axes = plt.subplots(self.ndim, figsize=(10, 3.5*self.ndim),
                                 sharex=True)
        samples = self.sampler.get_chain()
        labels = [r"$\beta$", r"$\epsilon$"]
        for i in range(self.ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], "r", alpha=0.3)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)

        axes[-1].set_xlabel("Step number")

        fout = self.handle_obs + '_emceeFlatChains.png'
        print('Saving chains: ' + fout)
        plt.savefig(fout)

        flat_samples = self.sampler.get_chain(discard=100, thin=15, flat=True)
        fig = corner.corner(flat_samples, labels=[r"$\beta$", r"$\epsilon$"],
                    show_titles=True, quantiles=[0.16, 0.84],
                    truths=[self.fs8/self.bs8, 1.0], truth_color='r')
        fout = self.handle_obs + '_emceeCorner.png'
        print('Saving corner: ' + fout)
        plt.savefig(fout)




    def log_likelihood(self, theta):
        fs8, bs8, sigmav, alpha, epsilon = theta
        alpha = 1.0
        G = 2 * beta / (3 + beta)
        xi0 = self.get_AP_multipole(self.xi0_APSpline, epsilon)
        xibar = self.get_AP_multipole(self.xibar_APSpline, epsilon)
        xi2 = self.get_AP_multipole(self.xi2_APSpline, epsilon)

        model = G * (xi0 - xibar)
        chi2 = np.dot(np.dot((xi2 - model), self.icov), xi2 - model)
        loglike = -1/2 * chi2 - np.log((2*np.pi)**(len(self.cov)/2)) * np.sum(np.linalg.eig(self.cov)[0])
        return loglike

    def log_prior(self, theta):
        beta, epsilon = theta

        if  0.1 < beta < 1.0 and 0.8 < epsilon < 1.2:
            return 0.0
        else:
            return -np.inf

    def log_probability(self, theta):
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(theta)


    def theory_multipoles(self, fs8, bs8, sigma_v, alpha_perp, alpha_par, s, mu):

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
            rescaled_r[i] = np.trapz((r[i] * alpha_par) * np.sqrt(1. + (1. - mus ** 2) *
                            (alpha_perp ** 2 / alpha_par ** 2 - 1)), mus)

        x = rescaled_r
        y1 = xi_r_func(r)
        y2 = self.delta_r(r)
        y3 = self.Delta_r(r)
        y4 = self.sv(r)

        # build rescaled interpolating functions using the relabelled separation vectors
        rescaled_xi_r = InterpolatedUnivariateSpline(x, y1, k=3)
        rescaled_delta_r = InterpolatedUnivariateSpline(x, y2, k=3, ext=3)
        rescaled_Delta_r = InterpolatedUnivariateSpline(x, y3, k=3, ext=3)
        rescaled_sv = InterpolatedUnivariateSpline(x, y4, k=3, ext=3)
        sigma_v = alpha_par * sigma_v
 

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
                sy = sigma_v * rescaled_sv_norm_func(rr) * self.iaH

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


            



