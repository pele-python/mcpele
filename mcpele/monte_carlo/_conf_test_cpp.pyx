# distutils: language = c++
# distutils: sources = conf_test.cpp


import sys

#===============================================================================
# Check spherical container
#===============================================================================

cdef class _Cdef_CheckSphericalContainer(_Cdef_ConfTest):
    """This class is the python interface for the c++ pele::CheckSphericalContainer configuration test class implementation
    """
    cdef cppCheckSphericalContainer* newptr
    def __cinit__(self, radius, ndim):
        self.thisptr = shared_ptr[cppConfTest](<cppConfTest*> new cppCheckSphericalContainer(radius, ndim))
        self.newptr = <cppCheckSphericalContainer*> self.thisptr.get()
        
class CheckSphericalContainer(_Cdef_CheckSphericalContainer):
    """This class is the python interface for the c++ CheckSphericalContainer implementation.
    """