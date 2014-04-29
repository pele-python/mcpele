# distutils: language = c++
cimport cython
import sys
from libcpp cimport bool as cbool
import numpy as np
cimport numpy as np
cimport pele.potentials._pele as _pele
cimport pele.optimize._pele_opt as _pele_opt 
from _pele_mc cimport cppConfTest,_Cdef_ConfTest

cdef extern from "mcpele/conf_test.h" namespace "mcpele":
    cdef cppclass cppCheckSphericalContainer "mcpele::CheckSphericalContainer":
        cppCheckSphericalContainer(double) except+

#===============================================================================
# Check spherical container
#===============================================================================

cdef class _Cdef_CheckSphericalContainer(_Cdef_ConfTest):
    """This class is the python interface for the c++ pele::CheckSphericalContainer configuration test class implementation
    """
    def __cinit__(self, radius):
        self.thisptr = <cppConfTest*>new cppCheckSphericalContainer(radius)
        
class CheckSphericalContainer(_Cdef_CheckSphericalContainer):
    """This class is the python interface for the c++ CheckSphericalContainer implementation.
    """