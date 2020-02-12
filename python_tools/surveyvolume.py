from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from python_tools.galaxycat import GalaxyCatalogue

fname = '/home/epaillas/data/eboss/v5/galaxy_cats/eBOSS_LRG_clustering_NGC_v5.dat.fits'

zmin=0.6
zmax=1.0
omega_m=0.307
h=0.6777

cat = GalaxyCatalogue(catalogue_file=fname, is_box=False, randoms=False,
boss_like=True, h=h, omega_m=omega_m, zmin=zmin, zmax=zmax, bin_write=False)

cat.getSurveyVolume()