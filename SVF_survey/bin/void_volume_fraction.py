import numpy as np
import healpy as hp
import sys
from astropy import units as u
from astropy.cosmology import Planck15, z_at_value, FlatLambdaCDM

if len(sys.argv) != 8:
    print('Some arguments are missing.')
    print('1) input_voids')
    print('2) random_catalogue')
    print('3) output_file')
    print('4) zmin')
    print('5) zmax')
    print('6) Om0')
    print('7) H0')
    sys.exit()

void_cat = sys.argv[1]
random_cat = sys.argv[2]
outfile = sys.argv[3]
zmin = float(sys.argv[4])
zmax = float(sys.argv[5])
Om0 = float(sys.argv[6])
H0 = float(sys.argv[7])

print('----------')
print('Running void_volume_fraction.py')
print('void_cat: ' + void_cat)
print('random_cat: ' + random_cat)
print('outfile: ' + outfile)
print('zmin: ' + repr(zmin))
print('zmax: ' + repr(zmax))
print('Om0: ' + repr(Om0))
print('H0: ' + repr(H0))

if '.npy' in void_cat:
    void_cat = np.load(void_cat)
else:
    void_cat = np.genfromtxt(void_cat)

if '.npy' in random_cat:
    random_cat = np.load(random_cat)
else:
    random_cat = np.genfromtxt(random_cat)

# Fiducial cosmology
cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Neff=Planck15.Neff, Tcmb0=Planck15.Tcmb0,
m_nu=Planck15.m_nu)

# Construct HEALPIX mask of randoms
nside = 512
npix = hp.nside2npix(nside)
mask = np.zeros(npix)
ind = hp.pixelfunc.ang2pix(nside, random_cat[:,1], random_cat[:,0], nest=False)
mask[ind] = 1

volfrac_list = []
npoints = 1000

for i in range(len(void_cat)):
#    print((repr(i/len(void_cat)*100))[:4] + '%', end='\r')
    rvoid = void_cat[i,3]
    xvoid = void_cat[i,0]
    yvoid = void_cat[i,1]
    zvoid = void_cat[i,2]

    ra = np.random.uniform(0, 2*np.pi, npoints)
    cosdec = np.random.uniform(-1, 1, npoints)
    dec = np.arccos(cosdec)
    dum = np.random.uniform(0, 1, npoints)
    r = rvoid * dum ** (1/3)

    x = r * np.sin(dec) * np.cos(ra) + xvoid
    y = r * np.sin(dec) * np.sin(ra) + yvoid
    z = r * np.cos(dec) + zvoid

    dis = np.sqrt(x**2 + y**2 + z**2)
    dec = np.arctan2(np.sqrt(x**2 + y**2), z)
    ra = np.arctan2(y, x)

    zhi = z_at_value(cosmo.comoving_distance, dis.min() * u.Mpc)
    zlo = z_at_value(cosmo.comoving_distance, dis.max() * u.Mpc)
    zgrid = np.logspace(np.log10(zhi), np.log10(zlo), 50)
    dgrid = cosmo.comoving_distance(zgrid)
    z = np.interp(dis, dgrid.value, zgrid)

    nin = 0
    for i in range(npoints):
        ind = hp.pixelfunc.ang2pix(nside, dec[i], ra[i], nest=False)
        if mask[ind] == 1 and zmin < z[i] < zmax:
            nin += 1

    volfrac = nin / npoints
    volfrac_list.append(volfrac)

if '.npy' in outfile:
    np.save(outfile, volfrac_list)
else:
    np.savetxt(outfile, volfrac_list)
