import numpy as np
import sys

if len(sys.argv) != 3:
    print('Some arguments are missing.')
    print('1) input_data')
    print('2) output_data')
    sys.exit()

print('----------')
print('Running angular_to_cartesian.py')

fin = sys.argv[1]
fout = sys.argv[2]

print('input_data: ' + fin)
print('output_data: ' + fout)

if '.npy' in fin:
    data = np.load(fin)
else:
    data = np.genfromtxt(fin)

ra = data[:,0]
dec = data[:,1]
dis = data[:,2]
z = data[:,3]

# Convert to comoving cartesian coordinates
xdata = np.array(dis * np.sin(dec) * np.cos(ra))
ydata = np.array(dis * np.sin(dec) * np.sin(ra))
zdata = np.array(dis * np.cos(dec))

xdata = np.reshape(xdata, (len(xdata), 1))
ydata = np.reshape(ydata, (len(ydata), 1))
zdata = np.reshape(zdata, (len(zdata), 1))
z = np.reshape(z, (len(z), 1))

data_out = np.hstack([xdata, ydata, zdata, z])

if '.npy' in fout:
    np.save(fout, data_out)
else:
    np.savetxt(fout, data_out)
