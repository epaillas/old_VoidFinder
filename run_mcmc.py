import numpy as np
import matplotlib.pyplot as plt
import python_tools.voidfitter as vf
from python_tools.utilities import FigureUtilities

handle_obs = '/Volumes/BlackIce/eboss/minerva_cats/void_stats/\
redshift_redshift/\
Galaxies_HOD_z0.57_Redshift_Redshift.SVF_recen_ovl0.5_MedianVoids.VG_CCF_rmu'

handle_mocks = '/Volumes/BlackIce/eboss/minerva_cats/void_stats/\
redshift_redshift/\
Galaxies_HOD_*_z0.57_Redshift_Redshift.SVF_recen_ovl0.5_MedianVoids.VG_CCF_rmu'

CaiModel = vf.CaiModel(handle_obs, handle_mocks)

#CaiModel.run_mcmc(niter=10000)
