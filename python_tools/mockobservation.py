
import numpy as np
import sys
import glob

handle_mocks = '/Volumes/BlackIce/eboss/minerva_cats/void_stats/\
real_real/\
Galaxies_HOD_*_z0.57_Real_Real.SVF_recen_ovl0.5_MedianVoids.VG_CCF_rmu'

handle_obs = '/Volumes/BlackIce/eboss/minerva_cats/void_stats/\
real_real/\
Galaxies_HOD_z0.57_Real_Real.SVF_recen_ovl0.5_MedianVoids.VG_CCF_rmu'

mock_files = sorted(glob.glob(handle_mocks))
xi_smu_list = []

for mock_file in mock_files:
    data = np.genfromtxt(mock_file)
    data[np.isnan(data)] = -1 # Set invalid values to -1
    data[data == np.inf] = -1 # Set invalid values to -1 
    xi_smu_list.append(data)

xi_smu_list = np.asarray(xi_smu_list)
xi_smu = np.mean(xi_smu_list, axis=0)

print(np.shape(xi_smu_list))
print(np.shape(xi_smu))

np.savetxt(handle_obs, xi_smu)

    







