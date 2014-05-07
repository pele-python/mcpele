# distutils: language = c++
cimport cython
import sys
from _pele_mc cimport cppTakeStep,_Cdef_TakeStep

cdef extern from "mcpele/takestep.h" namespace "mcpele":
    cdef cppclass cppRandomCoordsDisplacement "mcpele::RandomCoordsDisplacement":
        cppRandomCoordsDisplacement(size_t, size_t) except +
    cdef cppclass cppGaussianCoordsDisplacement "mcpele::GaussianCoordsDisplacement":
        cppGaussianCoordsDisplacement(size_t, size_t) except +

#===============================================================================
# RandomCoordsDisplacement
#===============================================================================

cdef class _Cdef_RandomCoordsDisplacement(_Cdef_TakeStep):
    """This class is the python interface for the c++ pele::RandomCoordsDisplacement take step class implementation
    """
    cdef cppRandomCoordsDisplacement* newptr
    def __cinit__(self, ndim, rseed):
        self.thisptr = <cppTakeStep*>new cppRandomCoordsDisplacement(ndim, rseed)
        self.newptr = <cppRandomCoordsDisplacement*> self.thisptr
        
class RandomCoordsDisplacement(_Cdef_RandomCoordsDisplacement):
    """This class is the python interface for the c++ RandomCoordsDisplacement implementation.
    """

#===============================================================================
# GaussianCoordsDisplacement
#===============================================================================

cdef class _Cdef_GaussianCoordsDisplacement(_Cdef_TakeStep):
    """This class is the python interface for the c++ pele::RandomCoordsDisplacement take step class implementation
    """
    cdef cppGaussianCoordsDisplacement* newptr
    def __cinit__(self, ndim, rseed):
        self.thisptr = <cppTakeStep*>new cppGaussianCoordsDisplacement(ndim, rseed)
        self.newptr = <cppGaussianCoordsDisplacement*> self.thisptr
    
class GaussianCoordsDisplacement(_Cdef_GaussianCoordsDisplacement):
    """This class is the python interface for the c++ RandomCoordsDisplacement implementation.
    """