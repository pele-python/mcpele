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
    def test_basic_harmonic(self):
        self.tsA = 
    
if __name__ == "__main__":
    unittest.main()
