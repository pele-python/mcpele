"""
This example computes the numerical value of Pi by Monte Carlo.
Points are sampled uniformly at random from an n-dimensional cube of
side length 2R. By measuring the fraction of points contained within
a ball of radius R contained in the cube, we determine Pi.
"""
from __future__ import division
import numpy as np
import copy
from scipy.special import gamma
from mcpele.monte_carlo import _BaseMCRunner
from mcpele.monte_carlo import NullPotential
from mcpele.monte_carlo import UniformRectangularSampling
from mcpele.monte_carlo import CheckSphericalContainerConfig

def volume_nball(radius, n):
    return np.power(np.pi, n / 2) * np.power(radius, n) / gamma(n / 2 + 1)

def get_pi(accepted_fraction, ndim):
    return np.power(2 ** ndim * accepted_fraction * gamma(ndim / 2 + 1), 2 / ndim)
    
class MC(_BaseMCRunner):
    def set_control(self, temp):
        self.set_temperature(temp)

class ComputePi(object):
    def __init__(self, ndim=2, nsamples=1e4):
        #
        self.ndim = ndim
        self.nsamples = nsamples
        #
        self.radius = 44
        self.potential = NullPotential()
        self.mc = MC(self.potential, np.ones(self.ndim), 1, self.nsamples)
        self.step = UniformRectangularSampling(42, self.radius)
        self.mc.set_takestep(self.step)
        self.mc.set_report_steps(0)
        self.conftest_check_spherical_container = CheckSphericalContainerConfig(self.radius)
        self.mc.add_conf_test(self.conftest_check_spherical_container)
        self.mc.set_print_progress()
        self.mc.run()
        self.p = self.mc.get_accepted_fraction()
        self.pi = get_pi(self.p, self.ndim)
        
if __name__ == "__main__":
    nsamples = 1e5
    ndim_ = []
    res = []
    for ndim in xrange(2, 16):
        print("computing pi in {} dimensions".format(ndim))
        c = ComputePi(ndim=ndim, nsamples=nsamples)
        res.append(c.pi)
        ndim_.append(ndim)
    for (i, p) in zip(ndim_, res):
        print("dimension", i)
        print("pi", p)
        print("pi / np.pi", p / np.pi)
