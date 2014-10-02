from _pele_mc cimport cppTakeStep,_Cdef_TakeStep, shared_ptr

cdef extern from "mcpele/random_coords_displacement.h" namespace "mcpele":
    cdef cppclass cppRandomCoordsDisplacement "mcpele::RandomCoordsDisplacement":
        cppRandomCoordsDisplacement(size_t, double) except +
        size_t get_seed() except +
        void set_generator_seed(const size_t) except +
        
cdef extern from "mcpele/random_coords_displacement_adaptive.h" namespace "mcpele":
    cdef cppclass cppRandomCoordsDisplacementAdaptive "mcpele::RandomCoordsDisplacementAdaptive":
        cppRandomCoordsDisplacementAdaptive(const size_t, const double) except +
        size_t get_seed() except +
        void set_generator_seed(const size_t) except +
        double get_stepsize() except +
    
cdef extern from "mcpele/gaussian_coords_displacement.h" namespace "mcpele":
    cdef cppclass cppGaussianCoordsDisplacement "mcpele::GaussianCoordsDisplacement":
        cppGaussianCoordsDisplacement(size_t, double) except +
        size_t get_seed() except +
        void set_generator_seed(const size_t) except +
        size_t get_count() except +
        
cdef extern from "mcpele/particle_pair_swap.h" namespace "mcpele":
    cdef cppclass cppParticlePairSwap "mcpele::ParticlePairSwap":
        cppParticlePairSwap(size_t, size_t) except +
        size_t get_seed() except +
        void set_generator_seed(const size_t) except +
        
cdef extern from "mcpele/adaptive_takestep.h" namespace "mcpele":
    cdef cppclass cppAdaptiveTakeStep "mcpele::AdaptiveTakeStep":
        cppAdaptiveTakeStep(shared_ptr[cppTakeStep], const size_t, const double,
                            const double, const double) except +
                            
cdef extern from "mcpele/take_step_pattern.h" namespace "mcpele":
    cdef cppclass cppTakeStepPattern "mcpele::TakeStepPattern":
        cppTakeStepPattern() except +
        void add_step(shared_ptr[cppTakeStep], const size_t)
        
cdef extern from "mcpele/take_step_probabilities.h" namespace "mcpele":
    cdef cppclass cppTakeStepProbabilities "mcpele::TakeStepProbabilities":
        cppTakeStepProbabilities(const size_t) except +
        void add_step(shared_ptr[cppTakeStep], const size_t)
