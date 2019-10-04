import numpy as np
import sys
from astropy.cosmology import Planck15, z_at_value, FlatLambdaCDM
import astropy.units as u

if len(sys.argv) != 5:
    print('Some arguments are missing.')
    print('1) input_data')
    print('2) output_data')
    print('3) Om0')
    print('4) H0')
    sys.exit()

print('----------')
print('Running cartesian_to_angular.py')

fin = sys.argv[1]
fout = sys.argv[2]
Om0 = float(sys.argv[3])
H0 = float(sys.argv[4])

print('input_data: ' + fin)
print('output_data: ' + fout)
print('Om0: ' + repr(Om0))
print('H0: ' + repr(H0))

if '.npy' in fin:
    data = np.load(fin)
else:
    data = np.genfromtxt(fin)

# Fiducial cosmology
cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Neff=Planck15.Neff, Tcmb0=Planck15.Tcmb0,
m_nu=Planck15.m_nu)

# Switch back to angular coordinates
dis = np.sqrt(data[:,0]**2 + data[:,1]**2 + data[:,2]**2)
dec = np.arctan2(np.sqrt(data[:,0]**2 + data[:,1]**2), data[:,2])
ra = np.arctan2(data[:,1], data[:,0])

rvoid = data[:,3]
ntrac = data[:,4]
nden = data[:,5]

zmin = z_at_value(cosmo.comoving_distance, dis.min() * u.Mpc)
zmax = z_at_value(cosmo.comoving_distance, dis.max() * u.Mpc)
zgrid = np.logspace(np.log10(zmin), np.log10(zmax), 50)
dgrid = cosmo.comoving_distance(zgrid)
z = np.interp(dis, dgrid.value, zgrid)

dis = np.reshape(dis, (len(dis), 1))
dec = np.reshape(dec, (len(dec), 1))
ra = np.reshape(ra, (len(ra), 1))
z = np.reshape(z, (len(z), 1))
rvoid = np.reshape(rvoid, (len(rvoid), 1))
ntrac = np.reshape(ntrac, (len(ntrac), 1))
nden = np.reshape(nden, (len(nden), 1))


data_out = np.hstack([ra, dec, z, rvoid, ntrac, nden])

if '.npy' in fout:
    np.save(fout, data_out)
else:
    np.savetxt(fout, data_out)
