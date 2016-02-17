from __future__ import division
import numpy as np
import pandas as pd

def write_2d_array_to_hf5(array, key, path):
    """
    this function can be used to dump any 2d array to a hf5 database
    use this to dump a trajectory to a database from python
    """
    assert array.ndim == 2
    nind , ncol = array.shape
    ind = [i for i in xrange(nind)]
    col = [i for i in xrange(ncol)]
    df = pd.DataFrame(np.array(array), index=ind, columns=col)
    df.to_hdf(path, key)

def read_hf5_to_2d_array(path, key):
    """
    read a 2d array from a hf5 database
    """
    df = pd.read_hdf(path, key)
    array = np.array(df.values)
    return array