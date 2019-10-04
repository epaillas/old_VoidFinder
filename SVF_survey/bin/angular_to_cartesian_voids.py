import numpy as np
import sys

if len(sys.argv) != 5:
    print('Some arguments are missing.')
    print('1) input_data')
    print('2) output_data')
    print('3) Om0')
    print('4) H0')
    sys.exit()

print('----------')
print('Running angular_to_cartesian.py')

fin = sys.argv[1]
fout = sys.argv[2]
Om = float(sys.argv[3])
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

ra = data[:,0]
dec = data[:,1]
z = data[:,2]
rvoid = data[:,3]
ntrac = data[:,4]
nden = data[:,5]

dis = np.array(cosmo.comoving_distance(z))

# Convert to comoving cartesian coordinates
xpos = np.array(dis * np.sin(dec) * np.cos(ra))
ypos = np.array(dis * np.sin(dec) * np.sin(ra))
zpos = np.array(dis * np.cos(dec))

xpos = np.reshape(xpos, (len(xpos), 1))
ypos = np.reshape(ypos, (len(ypos), 1))
zpos = np.reshape(zpos, (len(zpos), 1))
rvoid = np.reshape(rvoid, (len(rvoid), 1))
ntrac = np.reshape(ntrac, (len(ntrac), 1))
nden = np.reshape(nden, (len(nden), 1))

data_out = np.hstack([xpos, ypos, zpos, rvoid, ntrac, nden])

if '.npy' in fout:
    np.save(fout, data_out)
else:
    np.savetxt(fout, data_out)
