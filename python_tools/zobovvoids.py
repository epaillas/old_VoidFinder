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

class ZobovVoids:

    def __init__(self, handle, tracer_file, is_box=True,
                 box_size=1024):

        # file names
        self.handle = handle
        self.tracer_file = tracer_file
        self.tracer_unf = self.handle + '.dat.unf'
        self.is_box = is_box
        self.box_size = box_size

        self.ndiv = 2
        self.buffer_size = 0.1

        print('handle: ' + self.handle)
        print('tracer_file: ' + self.tracer_file)
        print('box_size: ' + str(self.box_size))
        print('is_box: ' + str(self.is_box))


        # convert input tracers to ZOBOV format
        self.ascii_to_bin()

        # run vozinit
        self.voz1b1()



    
    def ascii_to_bin(self):
        print('Converting tracer file to ZOBOV format...')

        if self.is_box:
            binpath = sys.path[0] + '/ZOBOV/bin/'
            cmd = [binpath + 'ascii_to_bin',
                   self.tracer_file,
                   self.tracer_unf
                  ]

            subprocess.call(cmd)

        return

    def vozinit(self):
        print('Running vozinit...')

        ndiv = 3
        buffer_size = 0.1

        if self.is_box:
            binpath = sys.path[0] + '/ZOBOV/bin/'
            cmd = [binpath + 'vozinit',
                   self.tracer_unf,
                   str(self.buffer_size),
                   str(self.box_size),
                   str(self.ndiv),
                   self.handle
                  ]
            subprocess.call(cmd)

        return

    def voz1b1(self):
        print('Running voz1b1...')

        if self.is_box:
            binpath = sys.path[0] + '/ZOBOV/bin/'

        for i in range(self.ndiv):
            for j in range(self.ndiv):
                for k in range(self.ndiv):
                    cmd = [binpath + 'voz1b1',
                           self.tracer_unf,
                           str(self.buffer_size),
                           str(self.box_size),
                           self.handle,
                           str(self.ndiv),
                           str(i),
                           str(j),
                           str(k)
                          ]

                    subprocess.call(cmd)

        return


