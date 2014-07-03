from _pele_mc cimport cppTakeStep,_Cdef_TakeStep, shared_ptr

cdef extern from "mcpele/takestep.h" namespace "mcpele":
    cdef cppclass cppRandomCoordsDisplacement "mcpele::RandomCoordsDisplacement":
        cppRandomCoordsDisplacement(size_t) except +
    cdef cppclass cppGaussianCoordsDisplacement "mcpele::GaussianCoordsDisplacement":
        cppGaussianCoordsDisplacement(size_t) except +