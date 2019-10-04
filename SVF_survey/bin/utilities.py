import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM, Planck15
import sys
from scipy.spatial import Delaunay
import healpy as hp


def fits_to_npy(fits_file):
    '''
    Read a BOSS/eBOSS fits galaxy catalogue
    that contains the RA, DEC and Z for
    each galaxy, and return it as a numpy array.
    '''
    hdul = fits.open(fits_file)
    data = hdul[1].data
    ra = data['RA']
    dec = data['DEC']
    z = data['Z']
    ra = (ra - 90)
    ra[ra > 360] -= 360
    ra[ra < 0] += 360
    ra = ra * (np.pi / 180)
    dec = (90 - dec) * (np.pi / 180)  
    ra = np.reshape(ra, (len(ra), 1))
    dec = np.reshape(dec, (len(dec), 1))
    z = np.reshape(z, (len(z), 1))
    
    return ra, dec, z

def z_to_dis(z, Om0, H0):
    '''
    Convert redshifts to distances using
    a fiducial cosmology.
    '''
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Neff=Planck15.Neff,
    Tcmb0=Planck15.Tcmb0, m_nu=Planck15.m_nu)
    dis = np.array(cosmo.comoving_distance(z))
    dis = np.reshape(dis, (len(dis), 1))
    return dis

def dis_to_z(dis, Om0, H0):
    '''
    Convert distances to redshifts using
    a fiducial cosmology.
    '''
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Neff=Planck15.Neff, Tcmb0=Planck15.Tcmb0,
    m_nu=Planck15.m_nu)
    zmin = z_at_value(cosmo.comoving_distance, dis.min() * u.Mpc)
    zmax = z_at_value(cosmo.comoving_distance, dis.max() * u.Mpc)
    zgrid = np.logspace(np.log10(zmin), np.log10(zmax), 50)
    dgrid = cosmo.comoving_distance(zgrid)
    z = np.interp(dis, dgrid.value, zgrid)
    return z

def sky_to_cartesian(ra, dec, z, Om0, H0):
    '''
    Convert sky to cartesian coordinates.
    '''
    dis = z_to_dis(z=z, Om0=Om0, H0=H0)
    x = dis * np.sin(dec) * np.cos(ra)
    y = dis * np.sin(dec) * np.sin(ra)
    z = dis * np.cos(dec)
    x = np.reshape(x, (len(x), 1))
    y = np.reshape(y, (len(y), 1))
    z = np.reshape(z, (len(z), 1))
    return x, y, z

def cartesian_to_sky(x, y, z, Om0, H0):
    '''
    Convert cartesian to sky coordinates
    '''
    dis = np.sqrt(data[:,0]**2 + data[:,1]**2 + data[:,2]**2)
    dec = np.arctan2(np.sqrt(data[:,0]**2 + data[:,1]**2), data[:,2])
    ra = np.arctan2(data[:,1], data[:,0])
    z = dis_to_z(dis=dis, Om0=Om0, H0=H0)
    ra = np.reshape(ra, (len(ra), 1))
    dec = np.reshape(dec, (len(dec), 1))
    z = np.reshape(z, (len(z), 1))
    return ra, dec, z
    

def delaunay_triangulation(points):
    '''
    Make a Delaunay triangulation out of 
    a set of points in cartesian coordinates.
    Returns the vertices of each tetraedron.
    '''
    triangulation = Delaunay(points)
    simplices = triangulation.simplices.copy()
    vertices = points[simplices]
    return vertices

def circumsphere(tetra):
    '''
    Calculate the circumcentre and circumradius
    for an input tetrahedron.
    '''
    x0, x1, x2, x3 = tetra
    A = []
    B = []
    A.append((x1 - x0).T)
    A.append((x2 - x0).T)
    A.append((x3 - x0).T)
    A = np.asarray(A)
    B = np.sum(A**2, axis=1)
    B = np.asarray(B)
    C = np.linalg.inv(A).dot(B)
    centre = x0 + 0.5 * C
    radius = 0.5 * np.sqrt(np.sum(C**2))
    return centre, radius

def construct_survey_mask(nside, ra, dec, zmin, zmax):
    '''
    Construct a HEALPix mask of the
    angular footprint of the survey.
    Also builds a mask of the pixels
    that lie at the survey boundaries.
    '''
    npix = hp.nside2npix(nside)
    mask = np.zeros(npix)
    ind = hp.pixelfunc.ang2pix(nside, dec, ra,
    nest=False)
    mask[ind] = 1
    indarray = [i for i in range(npix)]
    neigh = hp.pixelfunc.get_all_neighbours(nside, indarray, nest=False).T

    border = np.zeros(npix)
    for i in range(npix):
        if mask[i] == 0:
            count = 0
            for j in range(8):
                ind = neigh[i, j]
                if mask[ind] == 0:
                    count = count + 1
            if 0 < count <= 8:
                border[i] = 1

    return mask, border

def filter_by_mask(data, mask, nside, zlo, zhi):
    '''
    Uses an input HEALPix mask
    to filter data points that lie
    outside the mask. Points are also
    filtered by an input redshift range.
    '''
    ind = hp.pixelfunc.ang2pix(nside, data[:,1], data[:,0], nest=False)
    data = data[(mask[ind] == 1) & (zlo < data[:,2]) & (data[:,2] < zhi)]
    return data

def make_guard_particles(nden, ralo, rahi, declo, dechi,
                         zlo, zhi, Om0, H0, nside, mask, border):
    '''
    Constructs a set of guard particles
    at the survey boundaries
    '''
    dz = 5e-3
    ra, dec, z = make_random_sphere(nden, ralo, rahi, declo, dechi,
                                    zlo, zhi, Om0, H0)
    data = np.hstack([ra, dec, z])
    print(np.shape(data))
    print(np.shape(ra))
    print(np.shape(border))
    ind = hp.pixelfunc.ang2pix(nside, dec, ra, nest=False)
    print(np.shape(ind))
    print(np.shape(border[ind]))
    angCap = data[border[ind] == 1]
    redCap = data[mask[ind] == 1]

    angCap = [i for i in angCap if (zlo < i[3] < zhi)]
    redCap = [i for i in redCap if (zlo - dz < i[3] < zlo)
              or (zhi < i[3] < zhi + dz)]
    angCap = np.asarray(angCap)
    redCap = np.asarray(redCap)
    return angCap, redCap

def make_random_sphere(nden, ralo, rahi, declo, dechi,
                       zlo, zhi, Om0, H0):
    '''
    Constructs a sphere with a random
    distribution of particles with a 
    given number density (must be in (h/Mpc)^3).
    '''
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Neff=Planck15.Neff, Tcmb0=Planck15.Tcmb0,
    m_nu=Planck15.m_nu)
    dishi = cosmo.comoving_distance(zhi).value * cosmo.h
    dislo = cosmo.comoving_distance(zlo).value * cosmo.h
    vol = 4/3 * np.pi * (dishi**3)
    npoints = int(vol * nden)

    ralist = []
    declist = []
    zlist = []

    for i in range(npoints):
        ra = np.random.uniform(0, 2*np.pi)
        cosdec = np.random.uniform(-1, 1)
        dec = np.arccos(cosdec)
        u = np.random.uniform(0, 1)
        z = zhi * u ** (1/3)

        if (ralo < ra < rahi) and (declo < dec < dechi) and (zlo < z < zhi):
            ralist.append(ra)
            declist.append(dec)
            zlist.append(z)

    ralist = np.asarray(ralist).reshape(len(ralist), 1)
    declist = np.asarray(declist).reshape(len(declist), 1)
    zlist = np.asarray(zlist).reshape(len(zlist), 1)

    return ralist, declist, zlist

