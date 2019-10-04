# Sort void catalog by increasing order of radius

import numpy
import sys

# Input and output files
f = sys.argv[1]

print('-----------------------')
print('Running sort.py')
print('Input parameters:\n')
print('input_voids: ' + f)
print('output_voids: ' + f)

# Sort catalog by void radius (decreasing order)
catalog = numpy.genfromtxt(f)
catalog = catalog[numpy.argsort(catalog[:, -3])]
catalog = catalog[::-1]

print('\nNumber of voids: ' + repr(len(catalog)))
print('Rewriting void catalogue in decreasing order of void radius.')

# Write to file
fmt = '%10.3f %10.3f %10.3f %10i %10.3f'
numpy.savetxt(f, catalog, fmt=fmt)
