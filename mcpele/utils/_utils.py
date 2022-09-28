from __future__ import division
from builtins import range
import numpy as np
import pandas as pd

def write_2d_array_to_hdf5(array, key, path):
    """
    this function can be used to dump any 2d array to a hdf5 database
    use this to dump a trajectory to a database from python
    """
    assert array.ndim == 2
    nind , ncol = array.shape
    ind = [i for i in range(nind)]
    col = [i for i in range(ncol)]
    df = pd.DataFrame(np.array(array), index=ind, columns=col)
    df.to_hdf(path, key)

def read_hdf5_to_2d_array(path, key):
    """
    read a 2d array from a hdf5 database
    """
    df = pd.read_hdf(path, key)
    array = np.array(df.values)
    return array
