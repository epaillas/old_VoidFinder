
import numpy as np
import sys
import os
import glob
import subprocess
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



class Model1:
    '''
    Void-galaxy RSD model presented
    in Cai et al. (2016).
    '''

    def __init__(self, delta_r_file, xi_r_file, sv_file, xi_smu_file):

        self.delta_r_file = delta_r_file
        self.xi_r_file = xi_r_file
        self.sv_file = sv_file
        self.xi_smu_file = xi_smu_file


        print("Setting up void RSD model #1.")

        # cosmology for Minerva
        self.om_m = 0.285
        self.s8 = 0.828
        self.cosmo = Cosmology(om_m=self.om_m, s8=self.s8)

        self.eff_z = 0.57 # effective redshift for LRGs
        self.b = 2.3 # bias for LRGs

        self.growth = self.cosmo.get_growth(self.eff_z)
        self.f = self.cosmo.get_f(self.eff_z)
        self.s8norm = self.s8 * self.growth 

        eofz = np.sqrt((self.om_m * (1 + self.eff_z) ** 3 + 1 - self.om_m))
        self.iaH = (1 + self.eff_z) / (100. * eofz) 

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

        integral = np.zeros_like(self.r_for_xi)
        for i in range(len(integral)):
            integral[i] = quad(lambda x: self.xi_r(x) * x ** 2, 0, self.r_for_xi[i], full_output=1)[0]
        xibar_r = 3 * integral / self.r_for_xi ** 3
        self.xibar_r = InterpolatedUnivariateSpline(self.r_for_xi, xibar_r, k=1, ext=3)

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
        s, mu, xi_smu_obs = self.readCorrFile(self.xi_smu_file)
        s, self.xi0 = self._getMonopole(s, mu, xi_smu_obs)
        s, self.xi2 = self._getQuadrupole(s, mu, xi_smu_obs)

        monofunc = InterpolatedUnivariateSpline(s, self.xi0, k=3)
        integral = np.zeros_like(s)
        for i in range(len(integral)):
            integral[i] = quad(lambda x: monofunc(x) * x ** 2, 0, s[i], full_output=1)[0]
        self.xibar = 3 * integral / s ** 3

        self.s_for_xi = s
        self.mu_for_xi = mu

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


    def theory_multipoles(self, fs8, bs8, sigma_v, alpha_perp, alpha_para, s, mu):

        monopole = np.zeros(len(s))
        quadrupole = np.zeros(len(s))
        true_mu = np.zeros(len(mu))
        xi_model = np.zeros(len(mu))
        scaled_fs8 = fs8 / self.s8norm
        scaled_bs8 = bs8 / self.s8norm
        beta = scaled_fs8 / scaled_bs8
        print(beta)

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
        y5 = self.xibar_r(r)

        # build rescaled interpolating functions using the relabelled separation vectors
        rescaled_xi_r = InterpolatedUnivariateSpline(x, y1, k=3)
        rescaled_delta_r = InterpolatedUnivariateSpline(x, y2, k=3, ext=3)
        rescaled_Delta_r = InterpolatedUnivariateSpline(x, y3, k=3, ext=3)
        rescaled_sv = InterpolatedUnivariateSpline(x, y4, k=3, ext=3)
        sigma_v = alpha_para * sigma_v
        rescaled_xibar_r = InterpolatedUnivariateSpline(x, y5, k=3, ext=3)
 

        for i in range(len(s)):
            for j in range(len(mu)):
                true_sperp = s[i] * np.sqrt(1 - mu[j] ** 2) * alpha_perp
                true_spar = s[i] * mu[j] * alpha_para
                true_s = np.sqrt(true_spar ** 2. + true_sperp ** 2.)
                true_mu[j] = true_spar / true_s

                r = true_s * (1 + beta * rescaled_xibar_r(true_s) * true_mu[j]**2 / 3.)

                xi_model[j] = rescaled_xi_r(r) + beta/3 * rescaled_xibar_r(r) \
                              + beta * true_mu[j]**2 * (rescaled_xi_r(r) - rescaled_xibar_r(r))


            # build interpolating function for xi_smu at true_mu
            mufunc = InterpolatedUnivariateSpline(true_mu[np.argsort(true_mu)],
                                                  xi_model[np.argsort(true_mu)], k=3)
            
            # get multipoles
            xaxis = np.linspace(-1, 1, 1000)

            yaxis = mufunc(xaxis) / 2
            monopole[i] = simps(yaxis, xaxis)

            yaxis = mufunc(xaxis) * 5 / 2 * (3 * xaxis**2 - 1) / 2
            quadrupole[i] = simps(yaxis, xaxis)
            
            
        return monopole, quadrupole



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


            



