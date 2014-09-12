from __future__ import division
#import numpy as np

class CDFAccumulator(object):
    """
    Implemented with dict, could be slow, then use c++ version in
    cdf_accumulator.h
    """
    def __init__(self):
        self.total_number = 0
        self.data_x = dict()
    def add(self, inp):
        already_in = (inp in self.data_x)
        if already_in:
            self.data_x[inp] += 1
        else:
            self.data_x[inp] = 1
        self.total_number += 1
    def get_vecdata(self):
        x = self.data_x.keys()
        x = sorted(x)
        cdf_x = []
        remaining_x = self.total_number
        for xi in x:
            cdf_x.append(remaining_x / self.total_number)
            remaining_x -= self.data_x[xi]
            del self.data_x[xi]
        return x, cdf_x
    