import numpy as np
import matplotlib.pyplot as plt
from python_tools.caimodel import CaiModel
from python_tools.utilities import FigureUtilities

handle_obs = '/Volumes/BlackIce/eboss/minerva_cats/void_stats/\
real_redshift/\
Galaxies_HOD_z0.57_Real_Redshift.SVF_recen_ovl0.5_MedianVoids.VG_CCF_rmu'

handle_mocks = '/Volumes/BlackIce/eboss/minerva_cats/void_stats/\
real_redshift/\
Galaxies_HOD_*_z0.57_Real_Redshift.SVF_recen_ovl0.5_MedianVoids.VG_CCF_rmu'

model = CaiModel(handle_obs, handle_mocks)

model.run_mcmc(niter=5000, ncpu=3, nwalkers=32)
