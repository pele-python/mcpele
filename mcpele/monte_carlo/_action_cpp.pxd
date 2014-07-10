cimport pele.potentials._pele as _pele
from _pele_mc cimport cppAction,_Cdef_Action, shared_ptr
from libcpp cimport bool as cbool

cdef extern from "mcpele/actions.h" namespace "mcpele": 
    cdef cppclass cppAdjustStep "mcpele::AdjustStep":
        cppAdjustStep(double, double, size_t, size_t) except +
    cdef cppclass cppRecordEnergyHistogram "mcpele::RecordEnergyHistogram":
        cppRecordEnergyHistogram(double, double, double, size_t) except +
        _pele.Array[double] get_histogram() except +
        void print_terminal(size_t) except +
        double get_max() except +
        double get_min() except +
        double get_mean() except +
        double get_variance() except +
    cdef cppclass cppRecordEnergyTimeseries "mcpele::RecordEnergyTimeseries":
        cppRecordEnergyTimeseries(const size_t, const size_t) except +
        _pele.Array[double] get_time_series() except +
        void clear() except +
    cdef cppclass cppRecordLowestEValueTimeseries "mcpele::RecordLowestEValueTimeseries":
        cppRecordLowestEValueTimeseries(const size_t, const size_t,
            shared_ptr[_pele.cBasePotential], const size_t, _pele.Array[double]
            , const size_t) except +
        _pele.Array[double] get_time_series() except +
        void clear() except +