from libcpp cimport bool as cbool
from _pele_mc cimport cppConfTest, _Cdef_ConfTest, shared_ptr

cdef extern from "mcpele/check_spherical_container.h" namespace "mcpele":
    cdef cppclass cppCheckSphericalContainer "mcpele::CheckSphericalContainer":
        cppCheckSphericalContainer(double, size_t) except +
        
cdef extern from "mcpele/check_spherical_container_config.h" namespace "mcpele":
    cdef cppclass cppCheckSphericalContainerConfig "mcpele::CheckSphericalContainerConfig":
        cppCheckSphericalContainerConfig(double) except +

cdef extern from "mcpele/conf_test_OR.h" namespace "mcpele":
    cdef cppclass cppConfTestOR "mcpele::ConfTestOR":
        cppConfTestOR() except +
        void add_test(shared_ptr[cppConfTest]) except +
