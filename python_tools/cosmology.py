import numpy as np
from scipy.integrate import quad
from scipy.special import hyp2f1


class Cosmology:

    def __init__(self, om_m=0.308, h=0.676, s8=0.828):
        c = 299792.458
        om_l = 1.0 - om_m
        ztab = np.linspace(0, 4, 1000)
        rtab = np.zeros_like(ztab)
        for i in range(len(ztab)):
            rtab[i] = quad(lambda x: 0.01 * c / np.sqrt(om_m * (1 + x) ** 3 + om_l), 0, ztab[i])[0]

        self.h = h
        self.c = c
        self.om_m = om_m
        self.om_l = om_l
        self.ztab = ztab
        self.rtab = rtab
        self.s8 = s8

    # comoving distance in Mpc/h
    def get_comoving_distance(self, z):
        return np.interp(z, self.ztab, self.rtab)

    def get_redshift(self, r):
        return np.interp(r, self.rtab, self.ztab)

    def get_growth(self, eff_z):
        az = 1. / (1 + eff_z)
        growth = az ** 2.5 * np.sqrt(self.om_l + self.om_m * az ** (-3.)) * \
                hyp2f1(5. / 6, 3. / 2, 11. / 6, -(self.om_l * az ** 3.) / self.om_m) / \
                hyp2f1(5. / 6, 3. / 2, 11. / 6, -self.om_l / self.om_m)
        return growth

    def get_f(self, eff_z):
        f = ((self.om_m * (1 + eff_z)**3.) / (self.om_m * (1 + eff_z)**3 + self.om_l))**0.55
        return f

