cimport pele.potentials._pele as _pele
#from pele.potentials._pele cimport array_wrap_np
from _pele_mc cimport cppAction,_Cdef_Action, shared_ptr
from libcpp cimport bool as cbool
from libcpp.deque cimport deque

# cython has no support for integer template argument.  This is a hack to get around it
# https://groups.google.com/forum/#!topic/cython-users/xAZxdCFw6Xs
# Basically you fool cython into thinking INT2 is the type integer,
# but in the generated c++ code you use 2 instead.
# The cython code MyClass[INT2] will create c++ code MyClass<2>.
cdef extern from *:
    ctypedef int INT2 "2"    # a fake type
    ctypedef int INT3 "3"    # a fake type

cdef extern from "mcpele/record_energy_histogram.h" namespace "mcpele": 
    cdef cppclass cppRecordEnergyHistogram "mcpele::RecordEnergyHistogram":
        cppRecordEnergyHistogram(double, double, double, size_t) except +
        _pele.Array[double] get_histogram() except +
        void print_terminal() except +
        double get_max() except +
        double get_min() except +
        double get_mean() except +
        double get_variance() except +
        int get_count() except +

cdef extern from "mcpele/record_pair_dist_histogram.h" namespace "mcpele":
    cdef cppclass cppRecordPairDistHistogram "mcpele::RecordPairDistHistogram"[ndim]:
        cppRecordPairDistHistogram(_pele.Array[double], size_t, size_t, size_t) except +
        _pele.Array[double] get_hist_r() except +
        _pele.Array[double] get_hist_gr(double, size_t) except +
        size_t get_eqsteps() except +

cdef extern from "mcpele/record_scalar_timeseries.h" namespace "mcpele":
    cdef cppclass cppRecordScalarTimeseries "mcpele::RecordScalarTimeseries":
        _pele.Array[double] get_time_series() except +
        void clear() except +
        cbool moving_average_is_stable(size_t, double) except +

cdef extern from "mcpele/record_energy_timeseries.h" namespace "mcpele":
    cdef cppclass cppRecordEnergyTimeseries "mcpele::RecordEnergyTimeseries":
        cppRecordEnergyTimeseries(size_t, size_t) except +
        
cdef extern from "mcpele/record_lowest_evalue_timeseries.h" namespace "mcpele":
    cdef cppclass cppRecordLowestEValueTimeseries "mcpele::RecordLowestEValueTimeseries":
        cppRecordLowestEValueTimeseries(size_t, size_t,
            shared_ptr[_pele.cBasePotential], size_t, _pele.Array[double]
            , size_t) except +
    
cdef extern from "mcpele/record_displacement_per_particle_timeseries.h" namespace "mcpele":
    cdef cppclass cppRecordDisplacementPerParticleTimeseries "mcpele::RecordDisplacementPerParticleTimeseries":
        cppRecordDisplacementPerParticleTimeseries(size_t, size_t,
            _pele.Array[double], size_t) except +

cdef extern from "mcpele/record_vector_timeseries.h" namespace "mcpele":
    cdef cppclass cppRecordVectorTimeseries "mcpele::RecordVectorTimeseries":
        deque[_pele.Array[double]] get_time_series() except +
        void clear() except +
        size_t get_record_every() except +

cdef extern from "mcpele/record_coords_timeseries.h" namespace "mcpele":
    cdef cppclass cppRecordCoordsTimeseries "mcpele::RecordCoordsTimeseries":
        cppRecordCoordsTimeseries(size_t, size_t, size_t) except +
        _pele.Array[double] get_mean_coordinate_vector() except +
        _pele.Array[double] get_mean2_coordinate_vector() except +
        _pele.Array[double] get_variance_coordinate_vector() except +
        size_t get_count() except +
