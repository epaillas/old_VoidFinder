import numpy as np
import sys
from astropy.io import fits
from python_tools.cosmology import Cosmology
from scipy.io import FortranFile


class VoidCatalogue:

    def __init__(self, catalogue_file, is_box=True, box_size=1024.0,
                pos_cols=[0, 1, 2], omega_m=0.31, h=0.6777, verbose=True,
                bin_write=True, output_file=None):

        self.is_box = is_box

        if verbose:
            print('Loading void catalogue from file..')
            
        data = np.genfromtxt(catalogue_file)

        if self.is_box:
            self.box_size = box_size

            # for uniform box, position data is in Cartesian format
            self.x = data[:, pos_cols[0]]
            self.y = data[:, pos_cols[1]]
            self.z = data[:, pos_cols[2]]
            self.radius = data[:, 3]
            self.ntracers = data[:,4]
            self.nden = data[:,5]

        else:
            # position information is ra, dec and redshift
            self.ra = data[:, pos_cols[0]]
            self.dec = data[:, pos_cols[1]]
            self.redshift = data[:, pos_cols[2]]
            self.radius = data[:, 3]
            self.ntracers = data[:,4]
            self.nden = data[:,5]

            cosmo = Cosmology(omega_m=omega_m, h=h)
            self.dist = cosmo.get_comoving_distance(self.redshift)
            self.x = self.dist * np.cos(self.dec * np.pi / 180) * np.cos(self.ra * np.pi / 180)
            self.y = self.dist * np.cos(self.dec * np.pi / 180) * np.sin(self.ra * np.pi / 180)
            self.z = self.dist * np.sin(self.dec * np.pi / 180)

        self.x = self.x.reshape(len(self.x), 1)
        self.y = self.y.reshape(len(self.y), 1)
        self.z = self.z.reshape(len(self.z), 1)

        if bin_write:
            cout = np.c_[self.x, self.y, self.z, self.radius, self.ntracers,
                        self.nden]
            np.savetxt(output_file, cout)