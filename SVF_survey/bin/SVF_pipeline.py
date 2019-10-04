import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM, Planck15
import sys
from scipy.spatial import Delaunay
import healpy as hp
from utilities import *


# Import data and random fits catalogs
fname_data = '/home/epaillasv/data/boss/dr10/galaxy_cats/\
galaxy_DR10v8_CMASS_North.fits.gz'
fname_random = '/home/epaillasv/data/boss/dr10/galaxy_cats/\
galaxy_DR10v8_CMASS_North.fits.gz'
data_ra, data_dec, data_z = fits_to_npy(fname_data)
random_ra, random_dec, random_z = fits_to_npy(fname_random)

print('len(data): ' + repr(len(data_ra)))
print('len(random): ' + repr(len(random_ra)))

# Define a fiducial cosmology
H0 = 67.77
Om0 = 0.31

# Convert redshifts to distances using the
# fiducial cosmology
data_x, data_y, data_z = sky_to_cartesian(ra=data_ra, dec=data_dec,
                                          z=data_z, Om0=Om0, H0=H0)
random_x, random_y, random_z = sky_to_cartesian(ra=random_ra, dec=random_dec,
                                                z=random_z, Om0=Om0, H0=H0)

# Construct the survey mask
nside = 128
zmin = 0.43
zmax = 0.7
mask, border = construct_survey_mask(nside=nside, ra=random_ra,
dec=random_dec, zmin=zmin, zmax=zmax)

# Construct a random sphere with a number density
# equal to 5 times that of the data catalogue
nden = 5 * 3e-4 # In (h/Mpc)^3
ralo = random_ra.min()
rahi = random_ra.max()
declo = random_dec.min()
dechi = random_dec.max()
zhi = 0.73 # Overestimate to construct guard particles
zlo = 0.40
angCap, redCap = make_guard_particles(nden, ralo, rahi, declo, dechi,
                                      zlo, zhi, Om0, H0, nside, mask, border)


if plot:
    import matplotlib.pyplot as plt
    plt.style.use('enrique')

    # Plot footprint
    fig, axs = plt.subplots(1, 1, figsize=(7,5))
    axs.scatter(data_ra, data_dec, s=1.0)
    ax.set_xlabel('ra (deg)', fontsize=17)
    ax.set_ylabel('dec (deg)', fontsize=17)
    plt.tight_layout()
    plt.savefig('~/data/boss/footprint_sky.png', format='png')
    plt.show()

    # Plot cartesian projection
    fig, axs = plt.subplots(1, 3, figsize=(12,5))
    axs[0].scatter(data_x, data_y, s=1.0)
    axs[1].scatter(data_x, data_z, s=1.0)
    axs[2].scatter(data_y, data_z, s=1.0)
    ax[0].set_xlabel('x (Mpc)', fontsize=17)
    ax[0].set_ylabel('y (Mpc)', fontsize=17)
    plt.tight_layout()
    plt.savefig('~/data/boss/footprint_cartesian.png', format='png')



