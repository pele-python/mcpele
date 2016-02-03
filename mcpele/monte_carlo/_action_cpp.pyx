# distutils: language = c++
# distutils: sources = actions.cpp

cimport cython
cimport numpy as np

import numpy as np
import sys

#===============================================================================
# Record Energy Histogram
#===============================================================================        

cdef class _Cdef_RecordEnergyHistogram(_Cdef_Action):
    cdef cppRecordEnergyHistogram* newptr
    def __cinit__(self, min, max, bin, eqsteps):
        self.thisptr = shared_ptr[cppAction](<cppAction*> new cppRecordEnergyHistogram(min, max, bin, eqsteps))
        self.newptr = <cppRecordEnergyHistogram*> self.thisptr.get()
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def get_histogram(self):
        """returns the histogram array
        
        Returns
        -------
        numpy.array
            energy histogram
        """
        cdef _pele.Array[double] histi = self.newptr.get_histogram()
        cdef double *histdata = histi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] hist = np.zeros(histi.size())
        cdef size_t i
        for i in xrange(histi.size()):
            hist[i] = histdata[i]
              
        return hist
        
    def print_terminal(self):
        """draws histogram on the terminal
        """
        self.newptr.print_terminal()
    
    def get_bounds_val(self):
        """get energy boundaries of the histogram
        
        Returns
        -------
        Emax : double
            maximum energy of the histogram
        Emin : double
            minimum energy of the histogram
        """
        Emin = self.newptr.get_min()
        Emax = self.newptr.get_max()
        return Emin, Emax
    
    def get_mean_variance(self):
        """get mean and variance of the histogram
        
        Returns
        -------
        mean : double
            first moment of the distribution
        variance : double
            second central moment of the distribution
        """
        mean = self.newptr.get_mean()
        variance = self.newptr.get_variance()
        return mean, variance
    
    def get_count(self):
        """get number of entries in histogram
        
        Returns
        -------
        count : integer
            number of entries in histogram
        """
        count = self.newptr.get_count()
        return count
        
class RecordEnergyHistogram(_Cdef_RecordEnergyHistogram):
    """Bins energies into a resizable histogram
    
    This class is the Python interface for the c++ RecordEnergyHistogram implementation.
    
    .. warning :: :class:`RecordEnergyHistogram` should only start recording
                  entries when the system is equilibrated, set the number of steps
                  to skip with the ``eqsteps`` parameter.
    
    Parameters
    ----------
    min : double
        guess for the minimum energy expected
    max : double
        guess for the maximum energy expected
    bin : double
        choice for the bin size
    eqsteps: int
        number of iterations to skip before starting to record entries
        
    """

#===============================================================================
# Record Pair Distances Histogram
#===============================================================================
cdef class  _Cdef_RecordPairDistHistogram(_Cdef_Action):
    cdef cppRecordPairDistHistogram[INT2]* newptr2
    cdef cppRecordPairDistHistogram[INT3]* newptr3
    cdef cppRecordPairDistHistogramQuench[INT2]* newptr4
    cdef cppRecordPairDistHistogramQuench[INT3]* newptr5
    cdef _pele_opt.GradientOptimizer optimizer
    def __cinit__(self, boxvec, nr_bins, eqsteps, record_every, opt=None):
        ndim = len(boxvec)
        assert(ndim == 2 or ndim == 3)
        assert(len(boxvec)==ndim)
        cdef np.ndarray[double, ndim=1] bv
        if opt is None:
            self.quench = False
            if ndim == 2:
                bv = np.array(boxvec, dtype=float)
                self.thisptr = shared_ptr[cppAction](<cppAction*> new cppRecordPairDistHistogram[INT2](_pele.Array[double](<double*> bv.data, bv.size), nr_bins, eqsteps, record_every))
                self.newptr2 = <cppRecordPairDistHistogram[INT2]*> self.thisptr.get()
            else:
                bv = np.array(boxvec, dtype=float)
                self.thisptr = shared_ptr[cppAction](<cppAction*> new cppRecordPairDistHistogram[INT3](_pele.Array[double](<double*> bv.data, bv.size), nr_bins, eqsteps, record_every))
                self.newptr3 = <cppRecordPairDistHistogram[INT3]*> self.thisptr.get()
        else:
            self.quench = True
            self.optimizer = opt
            if ndim == 2:
                bv = np.array(boxvec, dtype=float)
                self.thisptr = shared_ptr[cppAction](<cppAction*> new cppRecordPairDistHistogramQuench[INT2](_pele.Array[double](<double*> bv.data, bv.size), nr_bins, eqsteps, record_every, self.optimizer.thisptr))
                self.newptr4 = <cppRecordPairDistHistogramQuench[INT2]*> self.thisptr.get()
            else:
                bv = np.array(boxvec, dtype=float)
                self.thisptr = shared_ptr[cppAction](<cppAction*> new cppRecordPairDistHistogramQuench[INT3](_pele.Array[double](<double*> bv.data, bv.size), nr_bins, eqsteps, record_every, self.optimizer.thisptr))
                self.newptr5 = <cppRecordPairDistHistogramQuench[INT3]*> self.thisptr.get()
        self.ndim = ndim
        
    def get_hist_r(self):
        """get array of :math:`r` values for :math:`g(r)` measurement
        
        Returns
        -------
        numpy.array
            array of :math:`r` values for :math:`g(r)` histogram
        """
        cdef _pele.Array[double] histi
        if self.ndim == 2 and not self.quench:
            histi = self.newptr2.get_hist_r()
        elif self.ndim == 3 and not self.quench:
            histi = self.newptr3.get_hist_r()
        elif self.ndim == 2 and self.quench:
            histi = self.newptr4.get_hist_r()
        elif self.ndim == 3 and self.quench:
            histi = self.newptr5.get_hist_r()
        cdef double *histdata = histi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] hist = np.zeros(histi.size())
        cdef size_t i
        for i in xrange(histi.size()):
            hist[i] = histdata[i]      
        return hist
    
    def get_hist_gr(self, number_density, nr_particles):
        """get array of :math:`g(r)` values for :math:`g(r)` measurement
        
        Returns
        -------
        numpy.array
            array of array of :math:`g(r)` values for :math:`g(r)`
        """
        cdef _pele.Array[double] histi
        if self.ndim == 2 and not self.quench:
            histi = self.newptr2.get_hist_gr(number_density, nr_particles)
        elif self.ndim == 3 and not self.quench:
            histi = self.newptr3.get_hist_gr(number_density, nr_particles)        
        elif self.ndim == 2 and self.quench:
            histi = self.newptr4.get_hist_gr(number_density, nr_particles) 
        elif self.ndim == 3 and self.quench:
            histi = self.newptr5.get_hist_gr(number_density, nr_particles) 
        cdef double *histdata = histi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] hist = np.zeros(histi.size())
        cdef size_t i
        for i in xrange(histi.size()):
            hist[i] = histdata[i]      
        return hist
    
    def get_eqsteps(self):
        """get number of equilibration steps
        
        Returns
        -------
        int
            number of equilibration steps
        """
        if self.ndim == 2 and not self.quench:
            return self.newptr2.get_eqsteps()
        elif self.ndim == 3 and not self.quench:
            return self.newptr3.get_eqsteps()
        elif self.ndim == 2 and self.quench:
            return self.newptr4.get_eqsteps()
        elif self.ndim == 3 and self.quench:
            return self.newptr5.get_eqsteps()
        else:
            raise Exception("_Cdef_RecordPairDistHistogram: boxdim fail")

class RecordPairDistHistogram(_Cdef_RecordPairDistHistogram):
    """Record a pair distribution function histogram
    
    This class is the Python interface for the c++ mcpele::RecordPairDistHistogram implementation.
    The pair correlation function (or radial distribution function) describes how the density of a
    system of particles varies as a function of distance from a reference particle. In simplest terms 
    it is a measure of the probability of finding a particle at a distance of :math:`r` away from a given 
    reference particle.
    
    Every time the action is called, it accumulates the present configuration into the same :math:`g(r)` histogram.
    The action function calls ``add_configuration`` which accumulates the current configuration into the :math:`g(r)` 
    histogram. The :math:`g(r)` histogram can be read out at any point after that.
     
    Parameters
    ----------
    boxvec : numpy.array
        array of box side lengths
    nr_bins : int
        number of bins for the :math:`g(r)` histogram
    eqsteps : int
        number of equilibration steps to be excluded from :math:`g(r)` computation
    record_every : int
        after ``eqsteps`` steps have been done, record every ``record_everyth`` steps
    opt : pele graident optimizer (optimal)
        If opt is passed, this will quench the snapshot of coords before
        accumulating the distances to the g(r) histogram. This is
        intended to give the quenched g(r) mentioned here:
        http://dx.doi.org/10.1063/1.449840
    """
    
#===============================================================================
# RecordEnergyTimeseries
#===============================================================================
        
cdef class _Cdef_RecordEnergyTimeseries(_Cdef_Action):
    cdef cppRecordScalarTimeseries* newptr
    def __cinit__(self, niter, record_every):
        self.thisptr = shared_ptr[cppAction](<cppAction*> new cppRecordEnergyTimeseries(niter, record_every))
        self.newptr = <cppRecordScalarTimeseries*> self.thisptr.get()
        
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def get_time_series(self):
        """get a energy time series array
        
        Returns
        -------
        numpy.array
            array containing the energy time series
        """
        cdef _pele.Array[double] seriesi = self.newptr.get_time_series()
        cdef double *seriesdata = seriesi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] series = np.zeros(seriesi.size())
        cdef size_t i
        for i in xrange(seriesi.size()):
            series[i] = seriesdata[i]
              
        return series
    
    def clear(self):
        """clear time series container
        
        deletes the entries in the c++ container
        """
        self.newptr.clear()
    
class RecordEnergyTimeseries(_Cdef_RecordEnergyTimeseries):
    """Record a time series of the energy
    
    This class is the Python interface for the c++ bv::RecordEnergyTimeseries 
    :class:`Action` class implementation.
    
    Parameters
    ----------
    niter: int, Deprecated
        expected number of steps (to preallocate)
    record_every : int
        interval every which the energy is recorded
    """
    
#===============================================================================
# RecordLowestEValueTimeseries
#===============================================================================
        
cdef class _Cdef_RecordLowestEValueTimeseries(_Cdef_Action):
    cdef cppRecordScalarTimeseries* newptr
    cdef ranvec
    def __cinit__(self, niter, record_every, _pele.BasePotential landscape_potential, boxdimension,
                  ranvec, lbfgsniter):
        cdef np.ndarray[double, ndim=1] ranvecc = ranvec
        self.thisptr = shared_ptr[cppAction](<cppAction*> new 
                 cppRecordLowestEValueTimeseries(niter, record_every,
                                                     landscape_potential.thisptr, boxdimension,
                                                     _pele.Array[double](<double*> ranvecc.data, ranvecc.size), lbfgsniter))
        self.newptr = <cppRecordScalarTimeseries*> self.thisptr.get()
        
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def get_time_series(self):
        """get the lowest eigenvalue time series array
        
        Returns
        -------
        numpy.array:
            array time series of the lowest eigenvalue
        """
        cdef _pele.Array[double] seriesi = self.newptr.get_time_series()
        cdef double *seriesdata = seriesi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] series = np.zeros(seriesi.size())
        cdef size_t i
        for i in xrange(seriesi.size()):
            series[i] = seriesdata[i]
              
        return series
    
    def clear(self):
        """clear time series container
        
        deletes the entries in the c++ container
        """
        self.newptr.clear()
    
class RecordLowestEValueTimeseries(_Cdef_RecordLowestEValueTimeseries):
    """Record lowest eigenvalue of the inherent structure
    
    This class is the Python interface for the c++ RecordLowestEValueTimeseries :class:`Action` 
    class implementation.
    The structure is quenched to a minimum energy configuration (its inherent structure) and
    the lowest eigenvalue is computed by the Rayleight-Ritz method for lowest eigenvalue (which
    computationally cheaper than the diagonalisation of the Hessian). The zero modes are
    orthogonalised through the Gram-Schmidt orthogonalisation procedure.
    
    Parameters
    ----------
    niter: int, Deprecated
        expected number of steps (to preallocate)
    record_every : int
        interval every which the energy is recorded
    landscape_potential : :class:`BasePotential <pele:pele.potentials.BasePotential>`
        potential associated with particles (so the underlying potential energy surface)
    boxdimension: int
        dimensionality of the space (dimensionality of box)
    ranvec : numpy.array
        random vector of length equal to the number of degrees of freedom [len(coords)],
        required by the Gram-Schmidt orthogonalisation procedure
    lbfgsniter : int
        maximum number of steps for the LBFG-S minimisation of the Rayleigh quotient
    """
    
#===============================================================================
# RecordDisplacementPerParticleTimeseries
#===============================================================================
        
cdef class _Cdef_RecordDisplacementPerParticleTimeseries(_Cdef_Action):
    cdef cppRecordScalarTimeseries* newptr
    cdef initial_coords
    def __cinit__(self, niter, record_every, initial_coords, boxdimension):
        cdef np.ndarray[double, ndim=1] initialc = initial_coords
        self.thisptr = shared_ptr[cppAction](<cppAction*> new 
                 cppRecordDisplacementPerParticleTimeseries(niter, record_every,
                                                            _pele.Array[double](<double*> initialc.data, initialc.size), 
                                                            boxdimension))
        self.newptr = <cppRecordScalarTimeseries*> self.thisptr.get()
        
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def get_time_series(self):
        """return a the root mean square displacement time series array
        
        Returns
        -------
        np.array
            root mean square displacement array
        """
        cdef _pele.Array[double] seriesi = self.newptr.get_time_series()
        cdef double *seriesdata = seriesi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] series = np.zeros(seriesi.size())
        cdef size_t i
        for i in xrange(seriesi.size()):
            series[i] = seriesdata[i]
              
        return series
    
    def clear(self):
        """clear time series container
        
        deletes the entries in the c++ container
        """
        self.newptr.clear()

class RecordDisplacementPerParticleTimeseries(_Cdef_RecordDisplacementPerParticleTimeseries):
    """Record time series of the average root mean square displacement per particle at each step
    
    This class is the Python interface for the c++ RecordDisplacementPerParticleTimeseries 
    :class:`Action` class implementation.
    
    Parameters
    ----------
    niter: int, Deprecated
        expected number of steps (to preallocate)
    record_every : int
        interval every which the energy is recorded
    initial_coords : numpy.array
        initial system coordinates, used to compute rms distance
    boxdimension: int
        dimensionality of the space (dimensionality of box)
    """
