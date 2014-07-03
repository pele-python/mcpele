from libcpp cimport bool as cbool
from _pele_mc cimport cppConfTest,_Cdef_ConfTest, shared_ptr

cdef extern from "mcpele/conf_test.h" namespace "mcpele":
    cdef cppclass cppCheckSphericalContainer "mcpele::CheckSphericalContainer":
        cppCheckSphericalContainer(double, size_t) except+