#!/bin/bash

#-------- Universal settings --------

ncores=3
model_number=2

root_real=/media/epaillasv/BlackIce/eboss/minerva_cats/void_stats/real_real/
root_redshift=/media/epaillasv/BlackIce/eboss/minerva_cats/void_stats/real_redshift/

xi_smu_obs="$root_redshift"Galaxies_HOD_z0.57_Real_Redshift.SVF_recen_ovl0.5_MedianVoids.VG_CCF_rmu
xi_smu_mocks="$root_redshift"Galaxies_HOD_*_z0.57_Real_Redshift.SVF_recen_ovl0.5_MedianVoids.VG_CCF_rmu
xi_r="$root_real"Galaxies_HOD_z0.57_Real_Real.SVF_recen_ovl0.5_MedianVoids.VG_CCF_monopole
delta_r="$root_real"Galaxies_HOD_z0.57_Real_Real.SVF_recen_ovl0.5_MedianVoids.VM_CCF_monopole
sv_r="$root_real"Galaxies_HOD_z0.57_Real_Real.SVF_recen_ovl0.5_MedianVoids.VG_CCF_losvel
covmat="$root_redshift"Galaxies_HOD_z0.57_Real_Redshift.SVF_recen_ovl0.5_MedianVoids_Model2_CovMat.npy

#-------- Run (do not modify below this point) --------
python $HOME/code/void_finder/run_mcmc.py \
--model_number "$model_number" \
--ncores "$ncores" \
--xi_smu_obs "$xi_smu_obs" \
--xi_smu_mocks "$xi_smu_mocks" \
--xi_r "$xi_r" \
--delta_r "$delta_r" \
--sv_r "$sv_r" \
--covmat "$covmat"

