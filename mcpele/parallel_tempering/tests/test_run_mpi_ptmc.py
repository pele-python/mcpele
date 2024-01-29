from __future__ import division
from builtins import range
from past.utils import old_div
import os
import tempfile
import numpy as np
import subprocess
import unittest
import logging


def read_Visits(fname):
    """HACKED"""
    data = np.genfromtxt(fname, delimiter="\t")
    return data[:, 0], data[:, 1]


class TestPTRun(unittest.TestCase):
    def test_heat_capacity(self):
        self.bdim = 3
        self.natoms = 4
        self.nprocs = 4
        # create a temporary directory using the context manager
        with tempfile.TemporaryDirectory() as tmpdir:
            subprocess.run(
                f"mpiexec --oversubscribe -n {self.nprocs} python "
                f"{os.path.dirname(os.path.abspath(__file__))}/_test_run_mpi_ptmc.py "
                f"{tmpdir}", shell=True, check=True)
            temperatures = np.genfromtxt(
                os.path.join(tmpdir, "temperatures"), delimiter="\t"
            )
            for i in range(self.nprocs):
                d = tmpdir + "/{}".format(i)
                pre = "Visits.his."
                files = os.listdir(d)
                ftlist = []
                for s in files:
                    if pre in s:
                        ftlist.append(s)
                timel = []
                for s in ftlist:
                    t = float(s[len(pre) :])
                    timel.append(t)
                max_t = np.amax(timel)

                ener, hist = read_Visits(d + "/" + pre + "{}".format(max_t))

                T = temperatures[i]

                average = np.average(ener, weights=hist)

                average2 = np.average(np.square(ener), weights=hist)

                cv = old_div((average2 - average**2), (T**2))
                cv_true = self.natoms * self.bdim / 2.0

                self.assertLess(
                    cv - cv_true,
                    0.1,
                    "failed for replica of rank {} cv = {}".format(i, cv),
                )


if __name__ == "__main__":
    logging.basicConfig(filename="ParallelTempering.log", level=logging.DEBUG)
    unittest.main()
