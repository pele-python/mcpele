from __future__ import division
import numpy as np
from mcpele.monte_carlo import _BaseMCRunner
from mcpele.monte_carlo import RandomCoordsDisplacement
from mcpele.monte_carlo import RecordPairDistHistogram

class MC(_BaseMCRunner):
    def set_control(self, temp):
        self.set_temperature(temp)

class ComputeGR():
    def __init__(self, boxdim=2, nr_particles=100, phi=0.4, nr_steps=1e6):
        self.boxdim = boxdim
        self.nr_particles = nr_particles
        self.phi = phi
        self.origin = np.zeros(self.nr_partilcles * self.boxdim)
        self.potential = Harmonic(self.origin, 42, bdim=self.boxdim, com=False)
        self.nr_steps = nr_steps
        self.mc = MC(self.potential, self.origin, 1, self.nr_steps)
        self.step = RandomCoordsDisplacement(42, 1, single=True)
        self.mc.set_takestep(self.step)
        #self.test = 
        #self.add_config_test(self.test)
        self.eqsteps = self.nr_steps / 2
        self.mc.set_report_steps(self.eq_steps)
        self.gr = RecordPairDistHistogram(self.boxvec, 500, self.eqsteps)
        self.mc.add_action(self.gr)
    def run(self):
        self.mc.set_print_progress()
        self.mc.run()