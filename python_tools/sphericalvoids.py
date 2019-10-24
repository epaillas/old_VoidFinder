import numpy as np
import sys
import os
import glob
import subprocess
from astropy.io import fits
from python_tools.cosmology import Cosmology
from python_tools.galaxycat import GalaxyCatalogue
from scipy.spatial import Delaunay
import healpy as hp
from scipy.io import FortranFile
import matplotlib.pyplot as plt

class SphericalVoids:

    def __init__(self, tracer_file, is_box=True, random_file='',
                 boss_like=False, pos_cols='0,1,2', box_size=1024.0,
                 omega_m=0.31, h=0.6777, mask_file='', zmin=0.43, zmax=0.7,
                 verbose=False, handle='', nside=512, delta_voids=0.2,
                 rvoidmax=100, ncores=1, steps='1,2,3,4'):

        steps = [int(i) for i in steps.split(',')]
        pos_cols = [int(i) for i in pos_cols.split(',')]

        # file names
        self.handle = handle
        self.tracer_file = tracer_file
        self.random_file = random_file
        self.tracer_unf = self.handle + '.dat.unf'
        self.random_unf = self.handle + '.ran.unf'
        self.centres_file = self.handle + '.cen.unf'
        self.voids_file = self.handle + '.SVF'
        self.recentred_file = self.voids_file + '_recen'
        self.mask_file = mask_file
        self.steps = steps
        self.pos_cols = pos_cols
        self.use_guards = True

        # void parameters
        self.delta_voids = delta_voids
        self.rvoidmax = rvoidmax

        # catalog data
        self.is_box = is_box
        self.box_size = box_size
        self.zmin = zmin
        self.zmax = zmax
        self.nside = nside

        # set cosmology
        self.omega_m = omega_m
        self.h = h
        self.cosmo = Cosmology(omega_m=omega_m)

        if 1 not in steps:
            if self.mask_file == '':
                sys.exit('Mask file not provided. Aborting...')
            else:
                self.mask = hp.read_map(self.mask_file, nest=False, verbose=False)
        else:
        # load tracers and find void centres
            self.tracers = GalaxyCatalogue(catalogue_file=tracer_file, is_box=is_box,
            box_size=box_size, randoms=False, boss_like=boss_like, omega_m=omega_m,
            h=h, bin_write=True, output_file=self.tracer_unf, pos_cols=pos_cols)
            
            if self.is_box == False:
                if random_file == '':
                    sys.exit('Random catalogue is missing. Aborting...')
                else:
                    self.randoms = GalaxyCatalogue(catalogue_file=random_file, is_box=self.is_box, 
                                                randoms=True, boss_like=boss_like, omega_m=omega_m,
                                                h=h, bin_write=True, output_file=self.random_unf,
                                                pos_cols=pos_cols)
                    
                if self.mask_file == '':
                    print('No mask file provided. Generating a rough mask...')
                    self.mask = self.make_survey_mask()
                else:
                    self.mask = hp.read_map(self.mask_file, nest=False, verbose=False)
                    
                self.mask_borders = self.get_mask_borders()

            self.get_circumcentres()

        # Grow spheres from the centres found in the previous step
        if 2 in steps:
            voids = self.grow_spheres(ncores=ncores)

        # Find better void centres by shifting the original positions
        if 3 in steps:
            voids = self.recentre_spheres(ncores=ncores)

        # Sort spheres in decreasing order of radius
        # and remove overlapping spheres
        if 4 in steps:
            voids = self.sort_spheres()

            # filter by void volume fraction
            if not self.is_box:
                voids = self.filter_by_volume_fraction(threshold=0.9)

            voids = self.overlap_filter(overlap=0.0)
            voids = self.overlap_filter(overlap=0.2)
            voids = self.overlap_filter(overlap=0.5)

        if (2 in steps) or (3 in steps) or (4 in steps):
            self.x = voids[:,0]
            self.y = voids[:,1]
            self.z = voids[:,2]
            self.radius = voids[:,3]
            self.ngal = voids[:,4]
            self.nden = voids[:,5]


    def concat_files(self, input_files, output_file):
        with open(output_file, 'w+') as outfile:
            for fname in input_files:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

    def in_ipynb():
        try:
            cfg = get_ipython().config 
            if cfg['IPKernelApp']['parent_appname'] == 'ipython-notebook':
                return True
            else:
                return False
        except NameError:
            return False


    def get_periodic_images(self, data):
        '''
        Find the relevant images of a 
        set of points in a box that
        has boundary conditions.
        '''
        images = []
        buffer = self.box_size / 10 # probably an overkill

        for point in data:
            condx = ((self.box_size - buffer) < point[0]) or (point[0] < buffer)
            condy = ((self.box_size - buffer) < point[1]) or (point[1] < buffer)
            condz = ((self.box_size - buffer) < point[2]) or (point[2] < buffer)

            if condx and condy and condz:
                shiftx = point[0] + np.copysign(self.box_size, (buffer - point[0]))
                shifty = point[1] + np.copysign(self.box_size, (buffer - point[1]))
                shiftz = point[2] + np.copysign(self.box_size, (buffer - point[2]))
                images.append([shiftx, shifty, shiftz])
            if condx:
                shiftx = point[0] + np.copysign(self.box_size, (buffer - point[0]))
                shifty = point[1]
                shiftz = point[2]
                images.append([shiftx, shifty, shiftz])
            if condy:
                shiftx = point[0]
                shifty = point[1] + np.copysign(self.box_size, (buffer - point[1]))
                shiftz = point[2]
                images.append([shiftx, shifty, shiftz])
            if condz:
                shiftx = point[0]
                shifty = point[1]
                shiftz = point[2] + np.copysign(self.box_size, (buffer - point[2]))
                images.append([shiftx, shifty, shiftz])

        images = np.asarray(images)
        return images
            

    def make_survey_mask(self):
        '''
        Constructs a binary HEALPix mask
        of the survey footprint.
        1 == inside survey
        0 == outside survey
        '''
        npix = hp.nside2npix(self.nside)
        mask = np.zeros(npix)
        theta = np.radians(90 - self.randoms.dec)
        phi = np.radians(self.randoms.ra)
        ind = hp.pixelfunc.ang2pix(self.nside, theta, phi, nest=False)
        mask[ind] = 1
        self.mask_file = self.handle + '_mask_nside{}.fits'.format(self.nside)
        hp.write_map(self.mask_file, mask, overwrite=True)
        return mask
    
    def get_mask_borders(self):
        '''
        Gets the boundary pixels from
        a HEALPix mask.
        1 == inside boundary pixel
        0 == outside boundary pixel
        '''
        npix = len(self.mask)
        nside = hp.npix2nside(npix)
        border = np.zeros(npix)
        ind = [i for i in range(npix)]
        neigh = hp.pixelfunc.get_all_neighbours(nside, ind, nest=False).T
                
        for i in range(npix):
            if self.mask[i] == 0:
                count = 0
                for j in range(8):
                    ind = neigh[i, j]
                    if self.mask[ind] == 1:
                        count = count + 1
                if 0 < count <= 8:
                    border[i] = 1
        return border

    def gen_random_sphere(self):
        dra = 2
        ddec = 2
        dz = 0.01
        nden = 5e-4

        ralo = self.randoms.ra.min() - dra
        rahi = self.randoms.ra.max() + dra
        declo = self.randoms.dec.min() - ddec
        dechi = self.randoms.dec.max() + ddec
        zhi = self.zmax + dz
        zlo = self.zmin - dz

        rlo = self.cosmo.get_comoving_distance(zlo)
        rhi = self.cosmo.get_comoving_distance(zhi)

        vol = 4/3 * np.pi * (rhi**3)
        npoints = int(vol * nden)

        ralist = []
        declist = []
        rlist = []
        zlist = []

        for i in range(npoints):
            ra = np.random.uniform(0, 2*np.pi)
            cosdec = np.random.uniform(-1, 1)
            dec = np.arccos(cosdec)
            u = np.random.uniform(0, 1)
            z = zhi * u ** (1/3)

            if (ralo < ra < rahi) and (declo < dec < dechi) and (zlo < z < zhi):
                r = cosmo.comoving_distance(z).value
                ralist.append(ra)
                declist.append(dec)
                rlist.append(r)
                zlist.append(z)

        ralist = np.asarray(ralist).reshape(len(ralist), 1)
        declist = np.asarray(declist).reshape(len(declist), 1)
        rlist = np.asarray(rlist).reshape(len(rlist), 1)
        zlist = np.asarray(zlist).reshape(len(zlist), 1)

        sphere = np.hstack([ralist, declist, rlist, zlist])

        return sphere

    def gen_guard_particles(self):
     
        sphere = self.gen_random_sphere()
        border_pix = self.get_mask_borders()

        ind = hp.pixelfunc.ang2pix(self.nside, sphere[:,1], sphere[:,0], nest=False)
        angCap = sphere[border_pix[ind] == 1]
        redCap = sphere[self.mask[ind] == 1]

        dz = 0.005
        angCap = [i for i in angCap if (self.zmin < i[3] < self.zmax)]
        redCap = [i for i in redCap if (self.zmin - dz < i[3] < self.zmin) or (self.zmax < i[3] < self.zmax + dz)]
        angCap = np.asarray(angCap).reshape((len(angCap), 3))
        redCap = np.asarray(redCap).reshape((len(redCap), 3))

        return angCap, redCap

    
    def get_circumcentres(self, radius_limit=300, bin_write=True):
        '''
        Find the centre of the circumspheres
        associated to an input catalogue of
        tetrahedra.
        '''

        print('Finding circumcentres of tetrahedra...')
        vertices = self.delaunay_triangulation()
        cenx, ceny, cenz, r = [], [], [], []
        
        for tetra in vertices:
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
            if radius < radius_limit:
                cenx.append(centre[0])
                ceny.append(centre[1])
                cenz.append(centre[2])
                r.append(radius)
                
        cenx = np.asarray(cenx)
        ceny = np.asarray(ceny)
        cenz = np.asarray(cenz)

        if not self.is_box:
            # Transform to sky coordinates
            cendis = np.sqrt(cenx ** 2 + ceny ** 2 + cenz ** 2)
            cenred = self.cosmo.get_redshift(cendis)
            cendec = np.arctan2(np.sqrt(cenx**2 + ceny**2), cenz)
            cenra = np.arctan2(ceny, cenx)

            # keep only those centres that are inside the survey
            ind = hp.pixelfunc.ang2pix(self.nside, cendec, cenra, nest=False)
            in_survey = (self.mask[ind] == 1) & (self.zmin < cenred) & (cenred < self.zmax)
            cenx = cenx[in_survey]
            ceny = ceny[in_survey]
            cenz = cenz[in_survey]

        else:
            # keep only those centres inside the box
            in_box = ((cenx >= 0) & (cenx <= self.box_size) &
                      (ceny >= 0) & (ceny <= self.box_size) &
                      (cenz >= 0) & (cenz <= self.box_size)
                     )
            cenx = cenx[in_box]
            ceny = ceny[in_box]
            cenz = cenz[in_box]

        cenx = cenx.reshape(len(cenx), 1)
        ceny = ceny.reshape(len(ceny), 1)
        cenz = cenz.reshape(len(cenz), 1)

        cout = np.hstack([cenx, ceny, cenz])

        print('{} centres found.'.format(len(cout)))

        if bin_write:
            f = FortranFile(self.centres_file, 'w')
            npoints = len(cout)
            f.write_record(npoints)
            f.write_record(cout)
            f.close()
        else:
            np.savetxt(self.centres_file, cout)

        self.centres = cout

        return
        

    def delaunay_triangulation(self, guards=False):
        '''eBOSS_LRG_NGC_v4.SVF.recen
        Make a Delaunay triangulation overeBOSS_LRG_NGC_v4.SVF.recen
        the cartesian positions of the tracers.eBOSS_LRG_NGC_v4.SVF.recen
        Returns the vertices of tetrahedra.
        '''
        x = self.tracers.x
        y = self.tracers.y
        z = self.tracers.z
        points = np.hstack([x, y, z])

        print(np.shape(points))
        
        if self.is_box == False and self.use_guards == True:
            angCap, redCap = self.gen_guard_particles()
            points = np.vstack([points, angCap, redCap])

        print(np.shape(points))
        
        # add periodic images if dealing with a box
        if self.is_box:
            images = self.get_periodic_images(points)
            points = np.vstack([points, images])

        triangulation = Delaunay(points)
        simplices = triangulation.simplices.copy()
        vertices = points[simplices]
        print('{} vertices found.'.format(len(vertices)))
        return vertices

    def grow_spheres(self, ncores=1):
        '''
        Grow spheres from an input 
        catalogue of centres.
        '''
        print('Proceeding to grow spheres...')
        if self.is_box:
            binpath = sys.path[0] + '/SVF_box/bin/'
            self.ngrid = 100
            cmd = ['mpirun', '-np', str(ncores), binpath + 'grow_spheres.exe',
                self.tracer_unf, self.centres_file, self.voids_file, str(self.box_size),
                str(self.delta_voids), str(self.rvoidmax), str(self.ngrid)]
        else:
            binpath = sys.path[0] + '/SVF_survey/bin/'
            self.gridmin = -5000
            self.gridmax = 5000
            cmd = ['mpirun', '-np', str(ncores), binpath + 'grow_spheres.exe',
                self.tracer_unf, self.random_unf, self.centres_file, self.voids_file,
                str(self.delta_voids), str(self.rvoidmax), str(self.gridmin),
                str(self.gridmax)]
        logfile = self.handle + '_grow_spheres.log'
        log = open(logfile, "w+")

        subprocess.call(cmd, stdout=log, stderr=log)

        if ncores > 1:
            files = glob.glob(self.voids_file + '.*')
            self.concat_files(input_files=files, output_file=self.voids_file)
            #subprocess.call(['rm'] + files)

        voids = np.genfromtxt(self.voids_file)
        return voids

    def recentre_spheres(self, ncores=1):
        '''
        Find better centres for an input
        void catalogue.
        '''
        print('Recentring spheres...')

        if self.is_box:
            binpath = sys.path[0] + '/SVF_box/bin/'
            self.ngrid = 100
            cmd = ['mpirun', '-np', str(ncores), binpath + 'recentring.exe',
                self.tracer_unf, self.voids_file, self.recentred_file,
                str(self.box_size), str(self.delta_voids), str(self.rvoidmax),
                str(self.ngrid)]
        else:
            binpath = sys.path[0] + '/SVF_survey/bin/'
            self.gridmin = -5000
            self.gridmax = 5000
            cmd = ['mpirun', '-np', str(ncores), binpath + 'recentring.exe',
                self.tracer_unf, self.random_unf, self.voids_file,
                self.recentred_file, str(self.delta_voids), str(self.rvoidmax),
                str(self.gridmin), str(self.gridmax)]
        logfile = self.handle + '_recentring.log'
        log = open(logfile, "w+")

        subprocess.call(cmd, stdout=log, stderr=log)

        if ncores > 1:
            files = glob.glob(self.recentred_file + '.*')
            self.concat_files(input_files=files, output_file=self.recentred_file)
            #subprocess.call(['rm'] + files)

        voids = np.genfromtxt(self.recentred_file)
        return voids

    def sort_spheres(self, fname='', radius_col=3):
        '''
        Sort an input void catalogue in
        decreasing order of radius.
        '''
        print('Sorting spheres by decreasing radius...')

        if fname == '':
            fname = self.recentred_file

        voids = np.genfromtxt(fname)
        voids = voids[np.argsort(voids[:, radius_col])]
        voids = voids[::-1]
        
        fmt = 4*'%10.3f ' + '%10i ' + '%10.3f '
        np.savetxt(fname, voids, fmt=fmt)

        return voids

    def overlap_filter(self, overlap=0.0):

        self.filtered_file = self.recentred_file + '_ovl{}'.format(str(overlap))

        if self.is_box:
            binpath = sys.path[0] + '/SVF_box/bin/'
            self.ngrid = 100
            cmd = [binpath + 'overlapping.exe', self.recentred_file, self.filtered_file,
                   str(self.box_size), str(overlap), str(self.ngrid)]
        else:
            binpath = sys.path[0] + '/SVF_survey/bin/'
            self.gridmin = -5000
            self.gridmax = 5000
            cmd = [binpath + 'overlapping.exe', self.recentred_file, self.filtered_file,
                   str(overlap), str(self.gridmin), str(self.gridmax)]

        logfile = self.handle + '_overlapping.log'
        log = open(logfile, "w+")

        subprocess.call(cmd, stdout=log, stderr=log)

        voids = np.genfromtxt(self.filtered_file)
        return voids


    def filter_by_volume_fraction(self, threshold=0.95):
        '''
        Filters voids by their volume fraction
        in the survey.
        '''
        print('Filtering voids by volume fraction...')
        voids = np.genfromtxt(self.recentred_file)
        volfrac = self.get_void_volume_fraction(fname=self.recentred_file)
        voids = voids[volfrac > threshold]

        self.recentred_file = self.recentred_file + '_vf{}'.format(str(threshold))
        fmt = 4*'%10.3f ' +  '%10i ' + '%10.3f '
        np.savetxt(self.recentred_file, voids, fmt=fmt)
        
        return voids


    def get_void_volume_fraction(self, fname='', pos_cols=[0,1,2],
                                 radius_col=3):
        '''
        Compute the fraction of the void
        that is contained within the survey
        footprint.
        '''
        if fname == '':
            fname = self.recentred_file

        voids = np.genfromtxt(fname)

        volfrac = []
        for void in voids:

            # Generate random points around each void
            npoints = 1000
            ra = np.random.uniform(0, 2*np.pi, npoints)
            cosdec = np.random.uniform(-1, 1, npoints)
            dec = np.arccos(cosdec)
            _ = np.random.uniform(0, 1, npoints)
            r = void[radius_col] * _ ** (1/3)

            # shift points to the global coordinate system
            x = r * np.sin(dec) * np.cos(ra) + void[pos_cols[0]]
            y = r * np.sin(dec) * np.sin(ra) + void[pos_cols[1]]
            z = r * np.cos(dec) + void[pos_cols[2]]

            # switch to sky coordinates
            dis = np.sqrt(x**2 + y**2 + z**2)
            dec = np.arctan2(np.sqrt(x**2 + y**2), z)
            ra = np.arctan2(y, x)
            redshift = self.cosmo.get_redshift(dis)

            # compute volume fraction
            nin = 0
            for i in range(npoints):
                ind = hp.pixelfunc.ang2pix(self.nside, dec[i], ra[i], nest=False)
                if self.mask[ind] == 1 and self.zmin < redshift[i] < self.zmax:
                    nin += 1

            volfrac.append(nin / npoints)

        volfrac = np.asarray(volfrac)
        return volfrac

    
