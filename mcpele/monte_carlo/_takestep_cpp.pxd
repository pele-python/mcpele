cimport pele.potentials._pele as _pele
from _pele_mc cimport cppTakeStep,_Cdef_TakeStep, shared_ptr

cdef extern from "mcpele/random_coords_displacement.h" namespace "mcpele":
    cdef cppclass cppRandomCoordsDisplacement "mcpele::RandomCoordsDisplacement":
        cppRandomCoordsDisplacement(size_t, double, double) except +
        size_t get_seed() except +
        void set_generator_seed(size_t) except +
        size_t get_count() except +
        void set_count(size_t) except +
        double get_stepsize() except +
        void set_stepsize(double) except +
    cdef cppclass cppRandomCoordsDisplacementAll "mcpele::RandomCoordsDisplacementAll":
        cppRandomCoordsDisplacementAll(size_t, size_t, size_t, double, double) except +
        size_t get_seed() except +
        void set_generator_seed(size_t) except +
        size_t get_count() except +
        void set_count(size_t) except +
        double get_stepsize() except +
        void set_stepsize(double) except +
    cdef cppclass cppRandomCoordsDisplacementSingle "mcpele::RandomCoordsDisplacementSingle":
        cppRandomCoordsDisplacementSingle(size_t, size_t, size_t, double, double) except +
        size_t get_seed() except +
        void set_generator_seed(size_t) except +
        size_t get_count() except +
        void set_count(size_t) except +
        double get_stepsize() except +
        void set_stepsize(double) except +

cdef extern from "mcpele/uniform_spherical_sampling.h" namespace "mcpele":
    cdef cppclass cppUniformSphericalSampling "mcpele::UniformSphericalSampling":
        cppUniformSphericalSampling(size_t, double) except +
        void set_generator_seed(size_t) except +
        void set_origin(_pele.Array[double]) except +

cdef extern from "mcpele/uniform_rectangular_sampling.h" namespace "mcpele":
    cdef cppclass cppUniformRectangularSampling "mcpele::UniformRectangularSampling":
        cppUniformRectangularSampling(size_t, _pele.Array[double]) except +
        void set_generator_seed(size_t) except +

cdef extern from "mcpele/gaussian_coords_displacement.h" namespace "mcpele":
    cdef cppclass cppGaussianTakeStep "mcpele::GaussianTakeStep":
        cppGaussianTakeStep(size_t, double, size_t) except +
        size_t get_seed() except +
        void set_generator_seed(size_t) except +
        size_t get_count() except +
        void set_count(size_t) except +
        double get_stepsize() except +
        void set_stepsize(double) except +
    cdef cppclass cppGaussianCoordsDisplacement "mcpele::GaussianCoordsDisplacement":
        cppGaussianCoordsDisplacement(size_t, double, size_t) except +
    cdef cppclass cppSampleGaussian "mcpele::SampleGaussian":
        cppSampleGaussian(size_t, double, _pele.Array[double]) except +

cdef extern from "mcpele/particle_pair_swap.h" namespace "mcpele":
    cdef cppclass cppParticlePairSwap "mcpele::ParticlePairSwap":
        cppParticlePairSwap(size_t, size_t, size_t) except +
        size_t get_seed() except +
        void set_generator_seed(size_t) except +

cdef extern from "mcpele/adaptive_takestep.h" namespace "mcpele":
    cdef cppclass cppAdaptiveTakeStep "mcpele::AdaptiveTakeStep":
        cppAdaptiveTakeStep(shared_ptr[cppTakeStep], size_t, double,
                            double, double) except +
        _pele.Array[size_t] get_counters() except +
        void set_counters(_pele.Array[size_t]) except +

cdef extern from "mcpele/take_step_pattern.h" namespace "mcpele":
    cdef cppclass cppTakeStepPattern "mcpele::TakeStepPattern":
        cppTakeStepPattern() except +
        void add_step(shared_ptr[cppTakeStep], size_t)

cdef extern from "mcpele/take_step_probabilities.h" namespace "mcpele":
    cdef cppclass cppTakeStepProbabilities "mcpele::TakeStepProbabilities":
        cppTakeStepProbabilities(size_t) except +
        void add_step(shared_ptr[cppTakeStep], double)
