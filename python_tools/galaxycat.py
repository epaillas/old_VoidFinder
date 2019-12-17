import numpy as np
import sys
from astropy.io import fits
from python_tools.cosmology import Cosmology
from scipy.io import FortranFile


class GalaxyCatalogue:

    def __init__(self, catalogue_file, is_box=True, box_size=1024.0, randoms=False, boss_like=False,
                pos_cols=[0, 1, 2], omega_m=0.31, h=0.6777, verbose=True, zmin=0, zmax=10,
                bin_write=True, output_file=None):

        self.is_box = is_box

        if not is_box and randoms and catalogue_file is None:
            sys.exit('Error: randoms file is missing.')

        if verbose:
            if randoms:
                print('Loading randoms data from file...')
            else:
                print('Loading galaxy data from file...')

        if boss_like:
            if self.is_box:
                sys.exit('Both boss_like and is_box cannot be simultaneously True. Aborting...')

            with fits.open(catalogue_file) as hdul:
                a = hdul[1].data
            for f in a.names:
                self.__dict__[f.lower()] = a.field(f)
  
            self.redshift = self.z

            # initialize Cartesian positions and observer distance
            cosmo = Cosmology(omega_m=omega_m, h=h)
            self.dist = cosmo.get_comoving_distance(self.redshift)
            self.x = self.dist * np.sin(self.dec * np.pi / 180) * np.cos(self.ra * np.pi / 180)
            self.y = self.dist * np.sin(self.dec * np.pi / 180) * np.sin(self.ra * np.pi / 180)
            self.z = self.dist * np.cos(self.dec * np.pi / 180)

        else:
            if '.npy' in catalogue_file:
                data = np.load(catalogue_file)
            else:
                try:
                    data = np.genfromtxt(catalogue_file)
                except:
                    sys.exit('Data format not recognized. Aborting...')

            if self.is_box:
                self.box_size = box_size

                # for uniform box, position data is in Cartesian format
                self.x = data[:, pos_cols[0]]
                self.y = data[:, pos_cols[1]]
                self.z = data[:, pos_cols[2]]

                if np.shape(data)[1] > 3:
                    print('Velocities included in data file...')
                    has_velocity = True
                    vel_cols=[pos_cols[0] + 3, pos_cols[1] + 3, pos_cols[2] + 3]
                    self.vx = data[:, vel_cols[0]]
                    self.vy = data[:, vel_cols[1]]
                    self.vz = data[:, vel_cols[2]]

            else:
                # position information is ra, dec and redshift
                self.ra = data[:, pos_cols[0]]
                self.dec = data[:, pos_cols[1]]
                self.redshift = data[:, pos_cols[2]]

                cosmo = Cosmology(omega_m=omega_m)
                self.dist = cosmo.get_comoving_distance(self.redshift)
                self.x = self.dist * np.sin(self.dec * np.pi / 180) * np.cos(self.ra * np.pi / 180)
                self.y = self.dist * np.sin(self.dec * np.pi / 180) * np.sin(self.ra * np.pi / 180)
                self.z = self.dist * np.cos(self.dec * np.pi / 180)

        self.x = self.x.reshape(len(self.x), 1)
        self.y = self.y.reshape(len(self.y), 1)
        self.z = self.z.reshape(len(self.z), 1)

        if has_velocity:
            self.vx = self.vx.reshape(len(self.vx), 1)
            self.vy = self.vy.reshape(len(self.vy), 1)
            self.vz = self.vz.reshape(len(self.vz), 1)

        # redshift cut
        if not self.is_box:
            ind = (self.redshift >= zmin) & (self.redshift <= zmax)
            self.x = self.x[ind]
            self.y = self.y[ind]
            self.z = self.z[ind]
            self.ra = self.ra[ind]
            self.dec = self.dec[ind]
            self.redshift = self.redshift[ind]

        if bin_write:
            if has_velocity:
                cout = np.hstack([self.x, self.y, self.z,
                                  self.vx, self.vy, self.vz])
            else:
                cout = np.hstack([self.x, self.y, self.z])

            f = FortranFile(output_file, 'w')
            nrows, ncols = np.shape(cout)
            f.write_record(nrows)
            f.write_record(ncols)
            f.write_record(cout)
            f.close()

class ProjectedGalaxyCatalogue:

    def __init__(self, catalogue_file, is_box=True, box_size=1024.0, randoms=False, boss_like=False,
                pos_cols=[0, 1], omega_m=0.31, h=0.6777, verbose=True, zmin=0, zmax=10,
                bin_write=True, output_file=None):

        self.is_box = is_box

        if not is_box and randoms and catalogue_file is None:
            sys.exit('Error: randoms file is missing.')

        if verbose:
            if randoms:
                print('Loading randoms data from file...')
            else:
                print('Loading galaxy data from file...')

        if boss_like:
            if self.is_box:
                sys.exit('Both boss_like and is_box cannot be simultaneously True. Aborting...')

            with fits.open(catalogue_file) as hdul:
                a = hdul[1].data
            for f in a.names:
                self.__dict__[f.lower()] = a.field(f)
  
            self.redshift = self.z

            # initialize Cartesian positions and observer distance
            cosmo = Cosmology(omega_m=omega_m, h=h)
            self.dist = cosmo.get_comoving_distance(self.redshift)
            self.x = self.dist * np.sin(self.dec * np.pi / 180) * np.cos(self.ra * np.pi / 180)
            self.y = self.dist * np.sin(self.dec * np.pi / 180) * np.sin(self.ra * np.pi / 180)

        else:
            if '.npy' in catalogue_file:
                data = np.load(catalogue_file)
            else:
                try:
                    data = np.genfromtxt(catalogue_file)
                except:
                    sys.exit('Data format not recognized. Aborting...')

            if self.is_box:
                self.box_size = box_size

                # for uniform box, position data is in Cartesian format
                self.x = data[:, pos_cols[0]]
                self.y = data[:, pos_cols[1]]

            else:
                # position information is ra, dec and redshift
                self.ra = data[:, pos_cols[0]]
                self.dec = data[:, pos_cols[1]]
                self.redshift = data[:, pos_cols[2]]

                cosmo = Cosmology(omega_m=omega_m)
                self.dist = cosmo.get_comoving_distance(self.redshift)
                self.x = self.dist * np.sin(self.dec * np.pi / 180) * np.cos(self.ra * np.pi / 180)
                self.y = self.dist * np.sin(self.dec * np.pi / 180) * np.sin(self.ra * np.pi / 180)

        self.x = self.x.reshape(len(self.x), 1)
        self.y = self.y.reshape(len(self.y), 1)

        # redshift cut
        if not self.is_box:
            ind = (self.redshift >= zmin) & (self.redshift <= zmax)
            self.x = self.x[ind]
            self.y = self.y[ind]
            self.z = self.z[ind]
            self.ra = self.ra[ind]
            self.dec = self.dec[ind]
            self.redshift = self.redshift[ind]

        if bin_write:
            cout = np.hstack([self.x, self.y])
            f = FortranFile(output_file, 'w')
            nrows, ncols = np.shape(cout)
            f.write_record(nrows)
            f.write_record(ncols)
            f.write_record(cout)
            f.close()