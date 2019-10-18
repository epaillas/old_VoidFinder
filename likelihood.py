from python_tools.mockstatistics import MockStatistics
import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner

def log_likelihood(theta, x, y, cov):
    beta = theta
    G = 2 * beta / (3 + beta)
    model = G * x
    # cov_dxi0, cov_xi2, cov_xi02, cov_xi20 = cov
    # cov_full = cov_xi2 + G**2 * cov_dxi0 - G*cov_xi02 - G*cov_xi20
    icov = np.linalg.inv(cov)
    chi2 = np.dot(np.dot((y - model), icov), y - model)
    loglike = -1/2 * chi2 - np.log((2*np.pi)**(len(cov)/2)) * np.sum(np.linalg.eig(cov)[0])
    return loglike

def log_prior(theta):
    beta = theta
    if 0.1 < beta < 0.7:
        return 0.0
    else:
        return -np.inf

def log_probability(theta, x, y, cov):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, cov)



handle = '/home/epaillasv/data/eboss/v5/2pcf/EZmock_eBOSS_LRG_SGC_v5_*.2pcr_rmu'
mocks = MockStatistics(handle)


xi0 = mocks.getMonopole(epsilon=1)
xi2 = mocks.getQuadrupole(epsilon=1)
xiav = mocks.getAverageMonopole(epsilon=1)
dxi0 = mocks.getDiffMonopole(epsilon=1)

cov_xi0 = mocks.getCovarianceMatrix(xi0)
cov_xi2 = mocks.getCovarianceMatrix(xi2)
cov_xiav = mocks.getCovarianceMatrix(xiav)
cov_dxi0 = mocks.getCovarianceMatrix(dxi0)
cov_xi02 = mocks.getCrossCovarianceMatrix(xi0, xi2)
cov_xi20 = mocks.getCrossCovarianceMatrix(xi2, xi0)
s = mock.bins

mean_xi0 = xi0.mean(axis=0)
mean_dxi0 = dxi0.mean(axis=0)
mean_xi2 = xi2.mean(axis=0)

# Visualize data for sanity check
fig = plt.figure()
ax = plt.subplot()

ax.plot(s, mean_xi0, 'mono')
ax.plot(s, mean_xi2, 'quad')

ax.legend()
plt.show()

# Cosmological parameters
Om0 = 0.31
Ode0 = 1 - Om0
bias_LRG = 2.30
eff_z_LRG = 0.703
f_LRG = ((Om0 * (1 + eff_z_LRG)**3.) / (Om0 * (1 + eff_z_LRG)**3 + Ode0))**0.55
fid_beta = f_LRG / bias_LRG
fid_G = 2 * fid_beta / (3 + fid_beta)
print('Fiducial beta: %0.3f' % fid_beta)

cov = cov_xi2 + fid_G**2 * cov_dxi0 - fid_G*cov_xi02 - fid_G*cov_xi20

nwalkers = 32
ndim = 1
beta_start = fid_beta
p0 = beta_start + 1e-4*np.random.randn(nwalkers, 1)

niter = 10000
data = (mean_dxi0, mean_xi2, cov)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=data)
sampler.run_mcmc(p0, niter, progress=True);


tau = sampler.get_autocorr_time()
flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)

fig = corner.corner(flat_samples, labels=[r'$f\sigma_8$', r'$b\sigma_8$'], show_titles=True)

# fig = plt.figure()
# ax = plt.subplot()

# bins = mocks.bins
# xi2 = mocks.quadrupole
# err = np.diag(mocks.quadrupole_cov)

# ax.errorbar(bins, np.mean(xi2, axis=0), np.sqrt(err))

# plt.show()

