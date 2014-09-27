from _pele_mc cimport cppTakeStep,_Cdef_TakeStep, shared_ptr

cdef extern from "mcpele/takestep.h" namespace "mcpele":
    cdef cppclass cppRandomCoordsDisplacement "mcpele::RandomCoordsDisplacement":
        cppRandomCoordsDisplacement(size_t, double) except +
        size_t get_seed() except +
        void set_generator_seed(const size_t) except +
    cdef cppclass cppRandomCoordsDisplacementAdaptive "mcpele::RandomCoordsDisplacementAdaptive":
        cppRandomCoordsDisplacementAdaptive(const size_t, const double, const double, const double, const double) except +
        size_t get_seed() except +
        void set_generator_seed(const size_t) except +
    cdef cppclass cppGaussianCoordsDisplacement "mcpele::GaussianCoordsDisplacement":
        cppGaussianCoordsDisplacement(size_t, double) except +
        size_t get_seed() except +
        void set_generator_seed(const size_t) except +
    cdef cppclass cppParticlePairSwap "mcpele::ParticlePairSwap":
        cppParticlePairSwap(size_t, size_t, size_t) except +
        size_t get_seed() except +
        void set_generator_seed(const size_t) except +
        
cdef extern from "mcpele/adaptive_takestep.h" namespace "mcpele":
    cdef cppclass cppAdaptiveTakeStep "mcpele::AdaptiveTakeStep":
        cppAdaptiveTakeStep(shared_ptr[cppTakeStep], const size_t) except +
