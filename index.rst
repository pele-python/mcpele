.. mcpele documentation master file, created by
   sphinx-quickstart on Sun Nov 30 12:26:21 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

mcpele : Monte Carlo Python Energy Landscape Explorer
+++++++++++++++++++++++++++++++++++++++++++++++++++++

Source code: https://github.com/pele-python/mcpele

Documentation: http://pele-python.github.io/mcpele/

Flexible and efficient Monte Carlo general purpose framework 
and MPI/mpi4py based Replica Exchange Method, built on the `pele <https://github.com/pele-python/pele>`_ 
foundations. mcpele provides a seamless integration of the
tools for energy landscape exploration built in pele. 
The package also acts as a plugin for the `Nested Sampling  <https://github.com/js850/nested_sampling>`_ project.

Through its c++ interface, mcpele makes Monte Carlo simulations available to 
researchers with little programming experience, without having to compromise
on efficiency. Furthermore mcpele abstracts each element of a Monte Carlo 
simulation eliminating the need for frequent code rewriting that experienced 
Monte Carlo developers typically go through, thus reducing the time required for
the implementation of an idea and reducing the occurrence of bugs.

.. figure:: diagram_fluid.png

  Figure 1: Diagramatic representation of the mcpele framework. On the right
  is a hard spheres fluid equilibrated by uniform sampling in a cubic box with
  periodic boundary conditions.

mcpele has been authored by Stefano Martiniani, Ken J Schrenk and Jacob Stevenson at the University of Cambridge.
The project is publicly available under the GNU general public licence.

Tutorials
-----------
.. toctree::
   :maxdepth: 3

   getting_started
   
Reference
---------

.. toctree::
   :maxdepth: 2
	
   BaseMCRunner
	
.. toctree::
   :maxdepth: 2
	
   TakeStep

.. toctree::
   :maxdepth: 2
   
   ConfTest

.. toctree::
   :maxdepth: 2
   
   AcceptTest

.. toctree::
   :maxdepth: 2
	
   Action

.. toctree::
   :maxdepth: 2
	
   Base_MPI_Parallel_Tempering

.. toctree::
   :maxdepth: 2
    
   MPI_Parallel_Tempering

Modules
+++++++
.. toctree::
   :maxdepth: 1

   monte_carlo
   parallel_tempering

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

