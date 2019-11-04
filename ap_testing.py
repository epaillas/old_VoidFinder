import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.integrate import quad

def gluedBins(bins, factor):
    '''
    Glues adjacent bins for an input
    array. The resulting array has 
    <factor> times fewer elements than
    the original one.
    '''

    glued_bins = []
    for i in range(0, len(bins), factor):
        glued_bins.append(np.mean(bins[i:i+factor]))
    return np.asarray(glued_bins)

def readCorrFile(fname):
    '''
    Function that reads the 3D correlation
    function data from an ascii file.
    Returns the correlation function in
    bins of s and mu, together with the
    corresponding s and mu vectors.
    '''
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

def getQuadrupole(xs, mu, xi_smu, epsilon=1):
    '''
    Computes the quadrupole from an input
    3D correlation function. Distorts
    input s and mu vectors according to the
    value of the epsilon AP parameter.
    '''
        
    qpara = 1 * epsilon ** (-2/3)
    qperp = epsilon * qpara

    true_sperp = s * np.sqrt(1 - mu**2) * qperp
    true_spara = s * mu * qpara
    true_s = np.sqrt(true_spara ** 2 + true_sperp ** 2)
    true_mu = true_spara / true_s

    quadr = np.zeros(xi_smu.shape[0])
    for j in range(xi_smu.shape[0]):
        mufunc = InterpolatedUnivariateSpline(true_mu, xi_smu[j, :], k=3)
        quadr[j] = quad(lambda x: mufunc(x) * 5 * (3. * x ** 2 - 1) / 2., 0, 1, full_output=1)[0]


    return true_s, quadr


fin = '/Users/epaillas/data/eboss/\
EZmock_eBOSS_LRG_NGC_v5_0001_Om307.SVF_recen_vf0.9_ovl0.5_sky_postOm307.VG_CCF_rmu'

# Read data
s, mu, xi_smu = readCorrFile(fin)

# Initialize axes
fig, ax = plt.subplots(1, figsize=(7,5))

# Plot quadrupole with 3 different values of AP parameter
epsilon_values = [0.95, 1.0, 1.05]

for epsilon in epsilon_values:
    # Compute quadrupole
    true_s, xi2 = getQuadrupole(s, mu, xi_smu, epsilon=epsilon)

    # Glue bins to reduce noise (default is 30 bins)
    true_s = gluedBins(true_s, factor=2)
    xi2 = gluedBins(xi2, factor=2)

    # Plot data
    ax.plot(true_s, xi2, label=r'$\epsilon = {}$'.format(epsilon))

plt.xlabel(r'$s\ [h^{-1}\rm{Mpc}]$', fontsize=14)
plt.ylabel(r'$\xi_2(s)$', fontsize=14)
ax.legend(fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.tick_params(axis='both', which='minor', labelsize=14)
plt.tight_layout()
plt.show()