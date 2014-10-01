from __future__ import division
import unittest
import logging
from mcpele.utils import CDFAccumulator
import numpy as np
from numpy.random import random_sample

class TestCDFAcc(unittest.TestCase):
    def setUp(self):
        self.nr_points = 100000
        self.low = - 42
        self.high = 42.42
    def test_basics(self):
        cdf = CDFAccumulator()
        for _ in xrange(self.nr_points):
            cdf.add(self.low + random_sample() * (self.high - self.low))
        x, cdf_x = cdf.get_vecdata()
        self.assertEqual(len(x), len(cdf_x))
        self.assertLessEqual(len(x), self.nr_points)
        self.assertLessEqual(x[-1], self.high)
        self.assertLessEqual(self.low, x[0])
        self.assertAlmostEqual(cdf_x[0], 1)
        self.assertAlmostEqual(cdf_x[-1], 1 / self.nr_points)
        for i in xrange(len(x)):
            self.assertAlmostEqual(cdf_x[i], (self.high - x[i]) / (self.high - self.low), delta = 1 / np.sqrt(self.nr_points))
        
if __name__ == "__main__":
    logging.basicConfig(filename = 'CDFAccumulator.log', level = logging.DEBUG)
    unittest.main()