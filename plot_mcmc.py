import numpy as np
import matplotlib.pyplot as plt
from python_tools.cosmology import Cosmology
import emcee
import corner

chain = '/Volumes/BlackIce/eboss/minerva_cats/\
void_stats/redshift_redshift/\
Galaxies_HOD_*_z0.57_Redshift_Redshift.SVF_recen_ovl0.5_MedianVoids.VG_CCF_rmu_emceeChain.h5'

from scipy.special import hyp2f1

# cosmology for Minerva
om_m = 0.285
s8 = 0.828
cosmo = Cosmology(om_m=om_m, s8=s8)

eff_z = 0.57 # effective redshift for LRGs
b = 2.3 # bias for LRGs
G = (2 * b) / (3 * b)

growth = cosmo.get_growth(eff_z)
f = cosmo.get_f(eff_z)
fs8 = f * s8 * growth
bs8 = b * s8 * growth

reader = emcee.backends.HDFBackend(chain)

# tau = reader.get_autocorr_time()
# burnin = int(2 * np.max(tau))
# thin = int(0.5 * np.min(tau))
# samples = reader.get_chain(discard=burnin, flat=True, thin=thin)
# log_prob_samples = reader.get_log_prob(discard=burnin, flat=True, thin=thin)
# log_prior_samples = reader.get_blobs(discard=burnin, flat=True, thin=thin)

# print("burn-in: {0}".format(burnin))
# print("thin: {0}".format(thin))
# print("flat chain shape: {0}".format(samples.shape))
# print("flat log prob shape: {0}".format(log_prob_samples.shape))
# print("flat log prior shape: {0}".format(log_prior_samples.shape))

# Corner
flat_samples = reader.get_chain(discard=50, thin=15, flat=True)
fig = corner.corner(flat_samples, labels=[r"$\beta$", r"$\epsilon$"],
            show_titles=True, quantiles=[0.16, 0.84],
            truths=[fs8/bs8, 1.0], truth_color='r')
#fout = self.handle_obs + '_emceeCorner.png'
#print('Saving corner: ' + fout)
plt.show()
#plt.savefig(fout)

# Quadrupole comparison: best-fit v/s true cosmology
fig, axs = plt.subplots(2, 1, figsize=(10, 5))

handle_obs = '/Volumes/BlackIce/eboss/minerva_cats/void_stats/\
redshift_redshift/\
Galaxies_HOD_*_z0.57_Redshift_Redshift.SVF_recen_ovl0.5_MedianVoids.VG_CCF_rmu'

handle_mocks = '/Volumes/BlackIce/eboss/minerva_cats/void_stats/\
redshift_redshift/\
Galaxies_HOD_*_z0.57_Redshift_Redshift.SVF_recen_ovl0.5_MedianVoids.VG_CCF_rmu'

CaiModel = vf.CaiModel(handle_obs, handle_mocks, mock_observation=True)

xi0_true, xibar_true, xi2_true = CaiModel.theory_multipoles()