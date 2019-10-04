import numpy as np
import sys

if len(sys.argv) != 5:
    print('Some arguments are missing.')
    print('1) input_file')
    print('2) output_file')
    print('3) volume_file')
    print('4) threshold')
    sys.exit()

infile = sys.argv[1]
outfile = sys.argv[2]
volfrac = sys.argv[3]
cut = float(sys.argv[4])

if '.npy' in infile:
    infile = np.load(infile)
else:
    infile = np.genfromtxt(infile)

if '.npy' in volfrac:
    volfrac = np.load(volfrac)
else:
    volfrac = np.genfromtxt(volfrac)

ind = volfrac >= cut

infile = infile[ind]

if '.npy' in outfile:
    np.save(outfile, infile)
else:
    np.savetxt(outfile, infile)
