#cython: boundscheck=False
#cython: wraparound=False

cimport pele.potentials._pele as _pele
from pele.potentials._pele cimport shared_ptr
from libcpp cimport bool as cbool

#===============================================================================
# mcpele::TakeStep
#===============================================================================

cdef extern from "mcpele/mc.h" namespace "mcpele":
    cdef cppclass cppTakeStep "mcpele::TakeStep"

cdef class _Cdef_TakeStep(object):
    """This class is the python interface for the c++ mcpele::TakeStep base class implementation
    """
    cdef shared_ptr[cppTakeStep] thisptr
    
#===============================================================================
# mcpele::AcceptTest
#===============================================================================

cdef extern from "mcpele/mc.h" namespace "mcpele":
    cdef cppclass cppAcceptTest "mcpele::AcceptTest"
        
cdef class _Cdef_AcceptTest(object):
    """This class is the python interface for the c++ mcpele::AcceptTest base class implementation
    """
    cdef shared_ptr[cppAcceptTest] thisptr
        
#===============================================================================
# mcpele::ConfTest
#===============================================================================

cdef extern from "mcpele/mc.h" namespace "mcpele":
    cdef cppclass cppConfTest "mcpele::ConfTest"

cdef class _Cdef_ConfTest(object):
    """This class is the python interface for the c++ mcpele::ConfTest base class implementation
    """
    cdef shared_ptr[cppConfTest] thisptr

#===============================================================================
# mcpele::Action
#===============================================================================

cdef extern from "mcpele/mc.h" namespace "mcpele":
    cdef cppclass cppAction "mcpele::Action"
        
cdef class _Cdef_Action(object):
    """This class is the python interface for the c++ mcpele::Action base class implementation
    """
    cdef shared_ptr[cppAction] thisptr

#===============================================================================
# mcpele::MC
#===============================================================================

cdef extern from "mcpele/mc.h" namespace "mcpele":
    cdef cppclass cppMC "mcpele::MC":
        cppMC(shared_ptr[_pele.cBasePotential], _pele.Array[double]&, double) except +
        void one_iteration() except +
        void run(size_t) except +
        void set_temperature(double) except +
        double get_temperature() except +
        void set_stepsize(double) except +
        void add_action(shared_ptr[cppAction]) except +
        void add_accept_test(shared_ptr[cppAcceptTest]) except +
        void add_conf_test(shared_ptr[cppConfTest]) except +
        void add_late_conf_test(shared_ptr[cppConfTest]) except +
        void set_takestep(shared_ptr[cppTakeStep]) except +
        void set_coordinates(_pele.Array[double]&, double) except +
        void reset_energy() except +
        double get_energy() except +
        _pele.Array[double] get_coords() except +
        _pele.Array[double] get_trial_coords() except +
        double get_accepted_fraction() except +
        size_t get_iterations_count() except +
        double get_conf_rejection_fraction() except +
        double get_E_rejection_fraction() except +
        size_t get_neval() except +
        double get_norm_coords() except +
        void set_report_steps(size_t) except +
        void abort() except +
        void set_print_progress() except +
        void enable_input_warnings() except+
        void disable_input_warnings() except +
        cbool get_success() except+

cdef class _Cdef_BaseMC(object):
    """This class is the python interface for the c++ mcpele::MC base class implementation
    """
    cdef shared_ptr[cppMC] thisptr
