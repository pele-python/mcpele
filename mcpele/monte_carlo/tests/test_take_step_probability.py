from __future__ import division
import numpy as np
from pele.potentials import Harmonic
from mcpele.monte_carlo import _BaseMCRunner, GaussianCoordsDisplacement
from mcpele.monte_carlo import TakeStepProbabilities, TakeStepPattern
from mcpele.monte_carlo import RandomCoordsDisplacement
import unittest

class MC(_BaseMCRunner):
    
    def set_control(self, temp):
        self.set_temperature(temp)
        
class TestTakeStepProbability(unittest.TestCase):
    
    def setUp(self):
        self.ndim = 42
        self.k = 100
        self.bdim = 2
        self.origin = np.zeros(self.ndim)
        self.potential = Harmonic(self.origin, self.k, bdim=self.bdim, com=False)
        self.potential_pattern = Harmonic(self.origin, self.k, bdim=self.bdim, com=False)
        self.temp = 1
        self.nr_steps = 1e4
        self.mc = MC(self.potential, self.origin, self.temp, self.nr_steps)
        self.mc_pattern = MC(self.potential_pattern, self.origin, self.temp, self.nr_steps)
        
    def test_frequencies(self):
        self.tsA = GaussianCoordsDisplacement(42, 1)
        self.tsA_pattern = GaussianCoordsDisplacement(42, 1)
        self.tsB = GaussianCoordsDisplacement(44, 2)
        self.tsB_pattern = GaussianCoordsDisplacement(44, 2)
        self.step = TakeStepProbabilities(42)
        self.step.add_step(self.tsA, 1)
        self.step.add_step(self.tsB, 3)
        self.step_pattern = TakeStepPattern()
        self.step_pattern.add_step(self.tsA_pattern, 1)
        self.step_pattern.add_step(self.tsB_pattern, 3)
        freqA = 1 / (1 + 3)
        freqB = 1 - freqA
        self.mc.set_takestep(self.step)
        self.mc.run()
        self.mc_pattern.set_takestep(self.step_pattern)
        self.mc_pattern.run()
        self.assertAlmostEqual(freqA, self.tsA.get_count() / self.nr_steps, delta=1e-2)
        self.assertAlmostEqual(freqB, self.tsB.get_count() / self.nr_steps, delta=1e-2)
        self.assertAlmostEqual(freqA, self.tsA_pattern.get_count() / self.nr_steps, delta=1e-2)
        self.assertAlmostEqual(freqB, self.tsB_pattern.get_count() / self.nr_steps, delta=1e-2)
        
class TestTakeStepProbabilityHarmoinc(unittest.TestCase):
    
    def setUp(self):
        self.box_dimension = 3
        self.nr_particles = 10
        self.k = 42
        self.nr_dof = self.box_dimension * self.nr_particles
        self.origin = np.zeros(self.nr_dof)
        self.potential = Harmonic(self.origin, self.k, bdim=self.box_dimension, com=True)
        self.temp = 1
        self.nr_steps = 1e4
        self.mc = MC(self.potential, self.origin, self.temp, self.nr_steps)
        self.take_step_A = RandomCoordsDisplacement(42, 1, single=True, nparticles=self.nr_particles, bdim=self.box_dimension)
        self.take_step_B = RandomCoordsDisplacement(44, 0.1, single=True, nparticles=self.nr_particles, bdim=self.box_dimension)
        self.step = TakeStepProbabilities(46)
        self.weight_A = 22
        self.weight_B = 78
        self.step.add_step(self.take_step_A, self.weight_A)
        self.step.add_step(self.take_step_B, self.weight_B)
        self.mc.set_takestep(self.step)
        self.frequency_step_A = self.weight_A / (self.weight_A + self.weight_B)
        self.frequency_step_B = self.weight_B / (self.weight_A + self.weight_B)
    
    def test_basic_harmonic(self):
        self.mc.run()
        self.assertAlmostEqual(self.frequency_step_A, self.take_step_A.get_count() / self.nr_steps, delta=1e-2)
        self.assertAlmostEqual(self.frequency_step_B, self.take_step_B.get_count() / self.nr_steps, delta=1e-2)
    
if __name__ == "__main__":
    unittest.main()
