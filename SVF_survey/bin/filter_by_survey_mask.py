"""
Filters an ra, dec input file by the
HEALPix map of the survey region.
"""
import numpy as np
import healpy as hp
import sys

if len(sys.argv) != 6:
    print('Some arguments are missing.')
    print('1) input_file')
    print('2) output_file')
    print('3) random_file')
    print('4) zmin')
    print('5) zmax')
    sys.exit()

print('----------')
print('Running filter_by_survey_mask.py')

infile = sys.argv[1]
outfile = sys.argv[2]
rand = sys.argv[3]
zmin = float(sys.argv[4])
zmax = float(sys.argv[5])

if '.npy' in infile:
    infile = np.load(infile)
else:
    infile = np.genfromtxt(infile)

if '.npy' in rand:
    rand = np.load(rand)
else:
    rand = np.genfromtxt(rand)

# Construct HEALPIX mask of dataoms
nside = 512 # Healpix resolution
npix = hp.nside2npix(nside) # Convert to number of pixels
mask = np.zeros(npix)
ind = hp.pixelfunc.ang2pix(nside, rand[:,1], rand[:,0], nest=False)
mask[ind] = 1

indarray = [i for i in range(npix)]
neigh = hp.pixelfunc.get_all_neighbours(nside, indarray, nest=False).T

border = np.zeros(npix)
for i in range(npix):
    if mask[i] == 0:
        count = 0
        for j in range(8):
            ind = neigh[i, j]
            if mask[ind] == 1:
                count = count + 1
        if 0 < count <= 8:
            border[i] = 1

# Keep data that is inside survey
ind = hp.pixelfunc.ang2pix(nside, infile[:,1], infile[:,0], nest=False)
infile = infile[mask[ind] == 1]
infile = np.asarray(infile)

# Keep data that is within the correct redshift range
infile = [i for i in infile if zmin < i[3] < zmax]
infile = np.asarray(infile)

if '.npy' in outfile:
    np.save(outfile, infile)
else:
    np.savetxt(outfile, infile)
