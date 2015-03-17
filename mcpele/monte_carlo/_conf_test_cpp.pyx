# distutils: language = c++
# distutils: sources = conf_test.cpp


import sys

#===============================================================================
# Check spherical container
#===============================================================================

cdef class _Cdef_CheckSphericalContainer(_Cdef_ConfTest):
    cdef cppCheckSphericalContainer* newptr
    def __cinit__(self, radius, ndim):
        self.thisptr = shared_ptr[cppConfTest](<cppConfTest*> new cppCheckSphericalContainer(radius, ndim))
        self.newptr = <cppCheckSphericalContainer*> self.thisptr.get()
        
class CheckSphericalContainer(_Cdef_CheckSphericalContainer):
    """Check that the system is within a spherical container
    
    This class is the Python interface for the c++ pele::CheckSphericalContainer 
    configuration test class implementation
    
    Parameters
    ----------
    radius : double
        radius of the spherical container, centered at **0**
    ndim : int
        dimensionality of the space (box dimensionality)
    """
    
#===============================================================================
# Check spherical container config
#===============================================================================

cdef class _Cdef_CheckSphericalContainerConfig(_Cdef_ConfTest):
    cdef cppCheckSphericalContainerConfig* newptr
    def __cinit__(self, radius):
        self.thisptr = shared_ptr[cppConfTest](<cppConfTest*> new cppCheckSphericalContainerConfig(radius))
        self.newptr = <cppCheckSphericalContainerConfig*> self.thisptr.get()
        
class CheckSphericalContainerConfig(_Cdef_CheckSphericalContainerConfig):
    """Check that the configuration point of the system is within a spherical container
    
    Parameters
    ----------
    radius : double
        radius of the spherical container, centered at **0**
    """
