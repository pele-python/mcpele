from __future__ import division
import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt
from pele.potentials import HS_WCA
from pele.optimize import LBFGS_CPP
from mcpele.monte_carlo import _BaseMCRunner
from mcpele.monte_carlo import RandomCoordsDisplacement
from mcpele.monte_carlo import RecordPairDistHistogram
from mcpele.monte_carlo import MetropolisTest

class MC(_BaseMCRunner):
    def set_control(self, temp):
        self.set_temperature(temp)

class ComputeGR():
    def __init__(self, boxdim=2, nr_particles=100, hard_phi=0.4,
                 nr_steps=1e6, epsilon=1, alpha=0.1, verbose=False):
        # Settings.
        np.random.seed(42)
        # Input parameters.
        self.boxdim = boxdim
        self.nr_particles = nr_particles
        self.hard_phi = hard_phi
        self.nr_steps = nr_steps
        self.epsilon = epsilon
        self.alpha = alpha
        self.verbose = verbose
        # Derived quantities.
        self.hard_radii = np.ones(self.nr_particles)
        def volume_nball(radius, n):
            return np.power(np.pi, n / 2) * np.power(radius, n) / gamma(n / 2 + 1)
        self.box_length = np.power(np.sum(np.asarray([volume_nball(r, self.boxdim) for r in self.hard_radii])) / self.hard_phi, 1 / self.boxdim)
        self.box_vector = np.ones(self.boxdim) * self.box_length
        # HS-WCA potential.
        self.potential = HS_WCA(use_periodic=True, use_cell_lists=True,
                                ndim=self.boxdim, eps=self.epsilon,
                                sca=self.alpha, radii=self.hard_radii,
                                boxvec=self.box_vector)
        # Initial configuration by minimization.
        self.nr_dof = self.boxdim * self.nr_particles
        self.x = np.random.uniform(-0.5 * self.box_length, 0.5 * self.box_length, self.nr_dof)
        optimizer = LBFGS_CPP(self.x, self.potential)
        optimizer.run()
        if not optimizer.get_result().success:
            print ("warning: minimization has not converged")
        self.x = optimizer.get_result().coords.copy()
        # Potential and MC rules.
        self.temperature = 1
        self.mc = MC(self.potential, self.x, self.temperature, self.nr_steps)
        self.step = RandomCoordsDisplacement(42, 1, single=True, nparticles=self.nr_particles, bdim=self.boxdim)
        if self.verbose:
            print ("initial MC stepsize")
            print self.step.get_stepsize()
        self.mc.set_takestep(self.step)
        self.eq_steps = self.nr_steps / 2
        self.mc.set_report_steps(self.eq_steps)
        self.gr_quench = RecordPairDistHistogram(self.box_vector, 50, self.eq_steps, self.nr_particles, optimizer=optimizer)
        self.gr = RecordPairDistHistogram(self.box_vector, 50, self.eq_steps, self.nr_particles)
        self.mc.add_action(self.gr_quench)
        self.mc.add_action(self.gr)
        self.test = MetropolisTest(44)
        self.mc.add_accept_test(self.test)
    def run(self):
        self.mc.set_print_progress()
        if not self.verbose:
            self.mc.disable_input_warnings()
        self.mc.run()
        if self.verbose:
            print ("adapted MC stepsize")
            print self.step.get_stepsize()
    def show_result(self):
        r = self.gr.get_hist_r()
        number_density = self.nr_particles / np.prod(self.box_vector)
        gr = self.gr.get_hist_gr(number_density, self.nr_particles)
        grq = self.gr_quench.get_hist_gr(number_density, self.nr_particles)
        plt.plot(r, gr, "o-", label="Equilibrium")
        plt.plot(r, grq, "x-", label="Quench")
        plt.xlabel(r"Distance $r$")
        plt.ylabel(r"Radial distr. function $g(r)$")
        plt.legend()
        plt.show()
        
if __name__ == "__main__":
    box_dimension = 2
    nr_particles = 100
    hard_volume_fraction = 0.4
    nr_steps = 1e5
    alpha = 0.48
    verbose = False
    simulation = ComputeGR(boxdim=box_dimension,
                           nr_particles=nr_particles,
                           hard_phi=hard_volume_fraction,
                           nr_steps=nr_steps,
                           alpha=alpha,
                           verbose=verbose)
    simulation.run()
    simulation.show_result()
