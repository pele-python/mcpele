import glob
import os
import sys
import subprocess

from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from numpy.distutils.misc_util import has_cxx_sources
import numpy as np
import pele

## Numpy header files 
numpy_lib = os.path.split(np.__file__)[0] 
numpy_include = os.path.join(numpy_lib, 'core/include') 

# find pele path
# note: this is used both for the c++ source files and for the cython pxd files,
# neither of which are "installed".  This should really point to the source directory.
# So this will only work if pele was built in-place
try:
    pelepath = os.path.dirname(pele.__file__)[:-5]
except:
    sys.stderr.write("WARNING: could't find path to pele\n")
    sys.exit()

def generate_cython():
    cwd = os.path.abspath(os.path.dirname(__file__))
    print("Cythonizing sources")
    p = subprocess.call([sys.executable,
                          os.path.join(cwd, 'cythonize.py'),
                          'mcpele', "-I %s/pele/potentials/" % pelepath],
                         cwd=cwd)
    if p != 0:
        raise RuntimeError("Running cythonize failed!")

generate_cython()

#
# compile fortran extension modules
#

class ModuleList:
    def __init__(self, **kwargs):
        self.module_list = []
        self.kwargs = kwargs
    def add_module(self, filename):
        modname = filename.replace("/", ".")
        modname, ext = os.path.splitext(modname)
        self.module_list.append(Extension(modname, [filename], **self.kwargs))

setup(name='mcpele', 
      version='0.1', 
      description="mcpele is a library of monte carlo and parallel tempering routines buil on the pele foundation",
      url='https://github.com/pele-python/mcpele',
      packages=["mcpele",
                "mcpele.monte_carlo",
                "mcpele.parallel_tempering",
                # add the test directories
                "mcpele.monte_carlo.tests",
                "mcpele.parallel_tempering.tests",
                ],
        )

#
# build the c++ files
#

include_sources_mcpele = ["source/" + f for f in os.listdir("source/") 
                   if f.endswith(".cpp")]
include_dirs = [numpy_include, "source"]

include_sources_pele = [pelepath+"/source/" + f for f in os.listdir(pelepath+"/source") 
                   if f.endswith(".cpp")]

depends_mcpele = [os.path.join("source/mcpele", f) for f in os.listdir("source/mcpele/") 
           if f.endswith(".cpp") or f.endswith(".h") or f.endswith(".hpp")]

depends_pele = [os.path.join(pelepath+"/source/pele", f) for f in os.listdir(pelepath+"/source/pele") 
                if f.endswith(".cpp") or f.endswith(".h") or f.endswith(".hpp")]

# note: on my computer (ubuntu 12.04 gcc version 4.6.3), when compiled with the
# flag -march=native I run into problems.  Everything seems to run ok, but when
# I run it through valgrind, valgrind complains about an unrecognized
# instruction.  I don't have a clue what is causing this, but it's probably
# better to be on the safe side and not use -march=native
#extra_compile_args = ['-I/home/sm958/Work/pele/source','-std=c++0x',"-Wall", "-Wextra", "-O3", '-funroll-loops']
# uncomment the next line to add extra optimization options

include_pele_source = '-I'+ pelepath + '/source'
extra_compile_args = [include_pele_source,'-std=c++0x',"-Wall", '-Wextra','-pedantic','-O3'] #,'-DDEBUG'

# note: to compile with debug on and to override extra_compile_args use, e.g.
# OPT="-g -O2 -march=native" python setup.py ...

include_sources_all = include_sources_mcpele + include_sources_pele

depends_all = depends_mcpele + depends_pele

cxx_modules = [
    Extension("mcpele.monte_carlo._pele_mc", 
              ["mcpele/monte_carlo/_pele_mc.cxx"] + include_sources_all,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends_all,
              ),
    Extension("mcpele.monte_carlo._monte_carlo_cpp", 
              ["mcpele/monte_carlo/_monte_carlo_cpp.cxx"] + include_sources_all,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends_all,
              ),
    Extension("mcpele.monte_carlo._takestep_cpp", 
              ["mcpele/monte_carlo/_takestep_cpp.cxx"] + include_sources_all,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends_all,
              ),
    Extension("mcpele.monte_carlo._accept_test_cpp", 
              ["mcpele/monte_carlo/_accept_test_cpp.cxx"] + include_sources_all,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends_all,
              ),
    Extension("mcpele.monte_carlo._conf_test_cpp", 
              ["mcpele/monte_carlo/_conf_test_cpp.cxx"] + include_sources_all,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends_all,
              ),
    Extension("mcpele.monte_carlo._action_cpp", 
              ["mcpele/monte_carlo/_action_cpp.cxx"] + include_sources_all,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends_all,
              ),
    Extension("mcpele.monte_carlo._nullpotential_cpp", 
              ["mcpele/monte_carlo/_nullpotential_cpp.cxx"] + include_sources_all,
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c++", depends=depends_all,
              ),
               ]
setup(ext_modules=cxx_modules,
      )
