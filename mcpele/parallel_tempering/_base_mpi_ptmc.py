from __future__ import division
from builtins import range
from builtins import object
import abc
import numpy as np
import os
from mpi4py import MPI
import copy
import logging
from future.utils import with_metaclass


class _MPI_Parallel_Tempering(with_metaclass(abc.ABCMeta, object)):
    """Abstract class for MPI Parallel Tempering calculations

    :class:`_MPI_Parallel_Tempering` implements all the basic MPI routines. The initialisation
    function and the method to find the swap pattern need to implemented in a specific child method.

    The replica exchange method (REM) or parallel tempering (PT) is a Monte Carlo scheme that targets
    the slow equilibration of systems characterised by large barriers in their free energy landscape.
    In REM :math:`n` replicas of the system are simulated simultaneously in the canonical (NVT) ensemble.
    These systems differ in temperature and, while the high temperature replicas explore the high energy
    regions of phase space, easily crossing large free energy barriers, the low temperature replicas
    explore the low lying regions of the energy landscape. In the hierarchical picture of the energy landscape,
    in which intra-funnel equilibration is fast and inter-funnel is slow, the high temperature replicas explore the
    energy landscape hopping between funnels, while the low temperature replicas explore individual funnels
    accurately. The idea of REM is to introduce moves that swap configurations between the different replicas,
    thus making the high energy regions available to the low temperature simulations and *vice versa*. It
    is not hard to see why this makes the exploration of configurational space more efficient, thus accelerating
    equilibration.

    The acceptance rule for parallel tempering is

    .. math:: P( x_i \Rightarrow x_j) = min \{ 1, \exp [- (\\beta_j - \\beta_i) (E_i - E_j)] \}

    The choice of temperature scale is very important for efficient equilibration to occur. From the acceptance rule
    we see that the acceptance will be suppressed exponentially by large differences in energy between the two
    configurations, just as for the Metropolis algorithm, but also by large differences in temperature. Thus, we
    must keep this product of differences as small as possible, to allow for efficient swapping. As a rule of thumb this
    can be achieved by using a scale of geometrically increasing temperatures.

    It is important to note that REM is a truly equilibrium Monte Carlo method: the microscopic equilibrium of each ensemble
    is not disturbed by the swaps, hence the ensemble averages for each replica are just as valid as for ordinary Monte Carlo
    simulations (unlike simulated annealing, for which ensemble averages are not well defined). In addition, ensemble averages
    for the extended ensemble can be obtained by histogram reweighting technique. We also highlight the fact that REM moves are
    very cheap to perform since they do not require additional energy evaluations.

    .. note :: An optimal Parallel Tempering strategy should make sure that all MCMC walks take roughly the same amount of time.
               Besides this fundamental consideration, note that root (rank=0) is not an evil master but rather an enlightened
               dictator that leads by example: root is responsible to assign jobs and control parameters (e.g. temperature) to
               the slaves but it also performs MCMC walks along with them. For this reason it might be optimal to give root a
               set of control parameters for which the simulation is leaner so that it can start doing its own things while the
               slaves finish their work.

    Parameters
    ----------
    mcrunner : :class:`_BaseMCrunner`
        object of :class:`_BaseMCrunner` that performs
        the MCMC walks
    Tmax : double
        maximum temperature to simulate (or equivalent control parameters)
    Tmin : double
        minimum temperature to simulate (or equivalent control parameters)
    max_ptiter : int
        maximum number of Parallel Tempering iterations
    pfreq : int
        frequency with which histogram and other information is dumped
        to a file
    skip : int
        number of parallel tempering iteration for which swaps should
        not be performed. Swaps should be avoided for instance while
        adjusting the step size
    print_status : bool
        choose whether to print MCrunner status at each iteration
    base_directory : string
        path to base directory where to save output

    Attributes
    ----------
    mcrunner : :class:`_BaseMCrunner`
        object of :class:`_BaseMCrunner` that performs
        the MCMC walks
    nprocs : int
        number of ranks
    rank : int
        MPI rank (identifier of the particular core), master has rank=0
    Tmax : double
        maximum temperature to simulate (or equivalent control parameters)
    Tmin : double
        minimum temperature to simulate (or equivalent control parameters)
    max_ptiter : int
        maximum number of Parallel Tempering iterations
    ptiter : int
        count of current parallel tempering iteration
    pfreq : int
        frequency with which histogram and other information is dumped
        to a file
    no_exchange_int : neg int
        this NEGATIVE number in :func:`exchange_pattern` means that no exchange should be attempted
    skip : int
        number of parallel tempering iteration for which swaps should
        not be performed. Swaps should be avoided for instance while
        adjusting the step size
    swap_accepted_count : int
        count of accepted swaps
    swap_rejected_count : int
        count of rejected swaps
    nodelist : list
        list of ranks (should be integers from 0 to ``nprocs``)
    base_directory : string
        path to base directory where to save output
    ex_outstream : stream
        stream to output exchanges informations
    initialised : bool
        records whether PT has been initialised
    print_status : bool
        choose whether to print MCrunner status at each iteration
    """

    def __init__(
        self,
        mcrunner,
        Tmax,
        Tmin,
        max_ptiter,
        pfreq=1,
        skip=0,
        print_status=True,
        base_directory=None,
    ):
        self.mcrunner = mcrunner
        self.comm = MPI.COMM_WORLD
        self.nprocs = self.comm.Get_size()  # total number of processes (replicas)
        self.rank = (
            self.comm.Get_rank()
        )  # this is the unique identifier for the process
        self.Tmax = Tmax
        self.Tmin = Tmin
        self.max_ptiter = max_ptiter
        self.ex_outstream = open("exchanges", "w")
        self.ptiter = 0
        self.print_status = print_status
        self.skip = (
            skip  # might want to skip the first few swaps to allow for equilibration
        )
        self.pfreq = pfreq
        self.no_exchange_int = (
            -12345
        )  # this NEGATIVE number in exchange pattern means that no exchange should be attempted
        self.initialised = False  # flag
        self.nodelist = [i for i in range(self.nprocs)]
        self.swap_accepted_count = 0
        self.swap_rejected_count = 0
        if base_directory is None:
            self.base_directory = os.path.join(os.getcwd(), "ptmc_results")
        else:
            self.base_directory = base_directory

        logging.info("Ready")

    @abc.abstractmethod
    def _find_exchange_buddy(self, Earray):
        """Abstract method to determines the exchange pattern, it needs to be overwritten.

        An exchange pattern array is constructed, filled with ``self.no_exchange_int`` which
        signifies that no exchange should be attempted. This value is replaced with the
        ``rank`` of the process with which to perform the swap if the swap attempt is successful.
        The exchange partner is then scattered to the other processes.

        Parameters
        ----------
        Earray : numpy.array
            array of energies (one from each core)
        """

    @abc.abstractmethod
    def _initialise(self):
        """Perform all the tasks required prior to starting the computation including initialising the output files"""

    @abc.abstractmethod
    def _print_data(self):
        """Abstract method responsible for printing and/or dumping the data, let it be printing the histograms or else"""

    @abc.abstractmethod
    def _print_status(self):
        """Abstract method responsible for printing and/or dumping the status, let it be printing the histograms or else"""

    @abc.abstractmethod
    def _print_exchanges(self):
        """Abstract method responsible for printing and/or dumping the number of exchanges per replica pair"""

    @abc.abstractmethod
    def _close_flush(self):
        """Abstract method responsible for printing and/or dumping the all streams at the end of the calculation"""

    def _test_convergence(self):
        """perform a convergence test, if it fails return new max_ptiter, else return self.max_ptiter"""
        return self.max_ptiter

    def one_iteration(self):
        """Perform one parallel tempering iteration

        Each PT iteration consists of the following steps:

        * set the coordinates
        * run the MCrunner for a predefined number of steps
        * collect the results (energy and new coordinates)
        * attempt an exchange
        """
        # set configuration and energy at which want to perform run
        self.mcrunner.set_config(np.array(self.config, dtype="d"), self.energy)
        # now run the MCMC walk
        self.mcrunner.run()
        # collect the results
        result = self.mcrunner.get_results()
        self.energy = result.energy
        self.config = np.array(result.coords, dtype="d")
        if self.ptiter >= self.skip:
            self._attempt_exchange()
            # print and increase parallel tempering count and test convergence
            if self.ptiter % self.pfreq == 0:
                self.max_ptiter = self._test_convergence()
                self._print_data()
            if self.print_status:
                self._print_status()
        self.ptiter += 1

    def run(self):
        """Run multiple single iterations, plus initialisation if MPI_PT has not been initialised yet"""
        if self.initialised is False:
            self._initialise()
        while self.ptiter < self.max_ptiter:
            if self.rank == self.nprocs - 1:
                logging.debug("Iteration {}".format(self.ptiter))
            self.one_iteration()
            if self.ptiter >= self.max_ptiter:
                self.max_ptiter = self._test_convergence()
                self._print_data()
        self._print_status()
        self._print_exchanges()
        self._close_flush()
        logging.info("Terminated")

    def _scatter_data(self, in_send_array, adim, dtype="d"):
        """Method to scatter data in equal ordered chunks among replica (it relies on the rank of the replica)

        In simple terms it requires that ``adim % nprocs = 0``.
        If root scatters an array of size ``nprocs`` then each core
        will receive the rank-th element of the array.

        Parameters
        ----------
        in_send_array : numpy.array
            incoming array (from the master) to be scattered
        adim : int
            size of the in_send_array
        dtype : dtype
            type of the elements of the array

        Returns
        -------
        recv_array : numpy.array
            array of length ``adim/nprocs``

        """
        if self.rank == 0:
            # process 0 is the root, it has data to scatter
            assert len(in_send_array) == adim
            assert adim % self.nprocs == 0
            send_array = np.array(in_send_array, dtype=dtype)
        else:
            # processes other than root do not send
            assert adim % self.nprocs == 0
            send_array = None

        recv_array = np.empty(int(adim / self.nprocs), dtype=dtype)
        self.comm.Scatter(send_array, recv_array, root=0)
        return recv_array

    def _scatter_single_value(self, send_array, dtype="d"):
        """Returns a single value from a scattered array for each replica (e.g. Temperature or swap partner)

        .. note :: send array must be of the same length as the number of processes

        Parameters
        ----------
        send_array : numpy.array
            incoming array (from the master) to be scattered
        dtype : dtype
            type of the elements of the array

        Returns
        -------
        dtype
            temperature or swap partner or else
        """
        if self.rank == 0:
            assert len(send_array) == self.nprocs

        T = self._scatter_data(send_array, self.nprocs, dtype=dtype)
        assert len(T) == 1
        return T[0]

    def _broadcast_data(self, in_data, adim, dtype="d"):
        """Identical arrays are broadcasted from root to all other processes

        Parameters
        ----------
        in_data : numpy.array
            incoming array (from the master) to be broadcasted
        adim : int
            size of the in_data array
        dtype : dtype
            type of the elements of the array

        Returns
        -------
        bcast_data : numpy.array
            array of length ``adim``

        """
        if self.rank == 0:
            bcast_data = np.array(in_data, dtype=dtype)
            assert len(bcast_data) == adim
        else:
            bcast_data = np.empty(adim, dtype=dtype)
        self.comm.Bcast(bcast_data, root=0)
        return bcast_data

    def _gather_data(self, in_send_array, dtype="d"):
        """Method to gather data in equal ordered chunks from replicas (it relies on the rank of the replica)

        .. note :: gather assumes that all the subprocess are sending the same amount of data to root, to send
                   variable amounts of data must use the MPI_gatherv directive

        Parameters
        ----------
        in_send_array : numpy.array
            incoming array (from each process) to be sent to master
        dtype : dtype
            type of the elements of the array

        Returns
        -------
        recv_array : numpy.array
            array of length ``len(in_send_array)) * nprocs``
        """
        in_send_array = np.array(in_send_array, dtype=dtype)
        if self.rank == 0:
            recv_array = np.empty(len(in_send_array) * self.nprocs, dtype=dtype)
        else:
            recv_array = None

        self.comm.Gather(in_send_array, recv_array, root=0)

        if self.rank != 0:
            assert recv_array is None

        return recv_array

    def _gather_energies(self, E):
        """gather energy of configurations from all processes

        Parameters
        ----------
        E : double
            energy of each process respectively

        Returns
        -------
        recv_Earray : numpy.array
            array of length ``nprocs``

        """
        send_Earray = np.array([E], dtype="d")
        recv_Earray = self._gather_data(send_Earray)
        return recv_Earray

    def _point_to_point_exchange_replace(self, dest, source, data):
        """swap data between two processes

        .. note :: the message sent buffer is replaced with the received message

        Parameters
        ----------
        dest : int
            rank of process with which to swap
        source : int
            rank of process with which to swap
        data : numpy.array
            array of data to exchage

        Returns
        -------
        data : numpy.array
            the send data buffer is replaced with the receive data
        """
        assert dest == source
        data = np.array(data, dtype="d")
        self.comm.Sendrecv_replace(data, dest=dest, source=source)
        return data

    def _exchange_pairs(self, exchange_buddy, data):
        """Return data from the pair exchange, otherwise return the data unaltered.

        .. warning :: the replica sends to exchange_partner and receives from it,
                      replacing source with self.rank would cause a deadlock

        Parameters
        ----------
        exchange_buddy : int
            rank of process with which to swap
        data : numpy.array
            array of data to exchage

        Returns
        -------
        data : numpy.array
            the send data buffer is replaced with the receive data
        """
        if exchange_buddy != self.no_exchange_int:
            # logging.debug("p-to-p exchange, old data {}".format(data))
            data = self._point_to_point_exchange_replace(
                exchange_buddy, exchange_buddy, data
            )
            # logging.debug("p-to-p exchange, new data {}".format(data))
            self.swap_accepted_count += 1
        else:
            self.swap_rejected_count += 1
        return data

    def _attempt_exchange(self):
        """This function brings together all the functions necessary to attempt a configuration swap

        This function is structures as follows:

        * root gathers the energies from the slaves
        * root decides who will swap with whom
        * root to each process the rank of its chosen partner (buddy)
        * processes exchange configuration and energy
        """
        # gather energies, only root will do so
        Earray = self._gather_energies(self.energy)
        if Earray is not None:
            logging.debug("Earray: {}".format(Earray))
        # find exchange pattern (list of exchange buddies)
        exchange_pattern = self._find_exchange_buddy(Earray)
        # now scatter the exchange pattern so that everybody knows who their buddy is
        exchange_buddy = self._scatter_single_value(
            np.array(exchange_pattern, dtype="d")
        )
        exchange_buddy = int(exchange_buddy)
        # attempt configurations swap
        self.config = self._exchange_pairs(
            exchange_buddy, np.array(self.config, dtype="d")
        )
        # swap energies
        E = self._exchange_pairs(exchange_buddy, np.array([self.energy], dtype="d"))
        assert len(E) == 1
        self.energy = E[0]
