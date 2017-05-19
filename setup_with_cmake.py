import os
import sys
import subprocess
import shutil
import argparse
import shlex

import numpy as np
from distutils import sysconfig
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from numpy.distutils.command.build_ext import build_ext as old_build_ext

import pele

# Numpy header files
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

# extract the -j flag and pass save it for running make on the CMake makefile
# extract -c flag to set compiler
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("-j", type=int, default=4)
parser.add_argument("-c", "--compiler", type=str, default=None)
jargs, remaining_args = parser.parse_known_args(sys.argv)

# record c compiler choice. use unix (gcc) by default
# Add it back into remaining_args so distutils can see it also
idcompiler = None
if not jargs.compiler or jargs.compiler in ("unix", "gnu", "gcc"):
    idcompiler = "unix"
    remaining_args += ["-c", idcompiler]
elif jargs.compiler in ("intel", "icc", "icpc"):
    idcompiler = "intel"
    remaining_args += ["-c", idcompiler]

# set the remaining args back as sys.argv
sys.argv = remaining_args
print jargs, remaining_args
if jargs.j is None:
    cmake_parallel_args = []
else:
    cmake_parallel_args = ["-j" + str(jargs.j)]

#extra compiler args
cmake_compiler_extra_args=["-std=c++0x","-Wall", "-Wextra", "-pedantic", "-O3"]
if idcompiler.lower() == 'unix':
    cmake_compiler_extra_args += ['-march=native', '-flto', '-fopenmp']
else:
    cmake_compiler_extra_args += ['-xavx', '-ipo', '-qopenmp']


#
# Make the git revision visible.  Most of this is copied from scipy
#
# Return the git revision as a string
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION

def write_version_py(filename='mcpele/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM SCIPY SETUP.PY
git_revision = '%(git_revision)s'
"""
    GIT_REVISION = git_version()

    a = open(filename, 'w')
    try:
        a.write(cnt % dict(git_revision=GIT_REVISION))
    finally:
        a.close()
write_version_py()
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

# note: on my computer (ubuntu 12.04 gcc version 4.6.3), when compiled with the
# flag -march=native I run into problems.  Everything seems to run ok, but when
# I run it through valgrind, valgrind complains about an unrecognized
# instruction.  I don't have a clue what is causing this, but it's probably
# better to be on the safe side and not use -march=native
#extra_compile_args = ['-I/home/sm958/Work/pele/source','-std=c++0x',"-Wall", "-Wextra", "-O3", '-funroll-loops', "-fopenmp", "-mavx"]
# uncomment the next line to add extra optimization options

# note: to compile with debug on and to override extra_compile_args use, e.g.
# OPT="-g -O2 -march=native" python setup.py ...

cmake_build_dir = "build/cmake"

cxx_files = ["mcpele/monte_carlo/_pele_mc.cxx",
             "mcpele/monte_carlo/_monte_carlo_cpp.cxx",
             "mcpele/monte_carlo/_takestep_cpp.cxx",
             "mcpele/monte_carlo/_accept_test_cpp.cxx",
             "mcpele/monte_carlo/_conf_test_cpp.cxx",
             "mcpele/monte_carlo/_action_cpp.cxx",
             "mcpele/monte_carlo/_nullpotential_cpp.cxx",
             ]

def get_ldflags(opt="--ldflags"):
    """return the ldflags.  This was taken directly from python-config"""
    getvar = sysconfig.get_config_var
    pyver = sysconfig.get_config_var('VERSION')
    libs = getvar('LIBS').split() + getvar('SYSLIBS').split()
    libs.append('-lpython'+pyver)
    # add the prefix/lib/pythonX.Y/config dir, but only if there is no
    # shared library in prefix/lib/.
    if opt == '--ldflags':
        if not getvar('Py_ENABLE_SHARED'):
            libs.insert(0, '-L' + getvar('LIBPL'))
        if not getvar('PYTHONFRAMEWORK'):
            libs.extend(getvar('LINKFORSHARED').split())
    return ' '.join(libs)

# create file CMakeLists.txt from CMakeLists.txt.in
with open("CMakeLists.txt.in", "r") as fin:
    cmake_txt = fin.read()
# We first tell cmake where the include directories are
cmake_txt = cmake_txt.replace("__PELE_INCLUDE__", pelepath + "/source")
# note: the code to find python_includes was taken from the python-config executable
python_includes = [sysconfig.get_python_inc(),
                   sysconfig.get_python_inc(plat_specific=True)]
cmake_txt = cmake_txt.replace("__PYTHON_INCLUDE__", " ".join(python_includes))
if isinstance(numpy_include, basestring):
    numpy_include = [numpy_include]
cmake_txt = cmake_txt.replace("__NUMPY_INCLUDE__", " ".join(numpy_include))
cmake_txt = cmake_txt.replace("__PYTHON_LDFLAGS__", get_ldflags())
cmake_txt = cmake_txt.replace("__COMPILER_EXTRA_ARGS__", '\"{}\"'.format(" ".join(cmake_compiler_extra_args)))
# Now we tell cmake which librarires to build
with open("CMakeLists.txt", "w") as fout:
    fout.write(cmake_txt)
    fout.write("\n")
    for fname in cxx_files:
        fout.write("make_cython_lib(${CMAKE_SOURCE_DIR}/%s)\n" % fname)

def set_compiler_env(compiler_id):
    """
    set environment variables for the C and C++ compiler:
    set CC and CXX paths to `which` output because cmake
    does not alway choose the right compiler
    """
    env = os.environ.copy()
    if compiler_id.lower() in ("unix"):
        env["CC"] = subprocess.check_output(["which", "gcc"]).rstrip('\n')
        env["CXX"] = subprocess.check_output(["which", "g++"]).rstrip('\n')
        env["LD"] = subprocess.check_output(["which", "ld"]).rstrip('\n')
        env["AR"] = subprocess.check_output(["which", "ar"]).rstrip('\n')
    elif compiler_id.lower() in ("intel"):
        env["CC"] = subprocess.check_output(["which", "icc"]).rstrip('\n')
        env["CXX"] = subprocess.check_output(["which", "icpc"]).rstrip('\n')
        env["LD"] = subprocess.check_output(["which", "xild"]).rstrip('\n')
        env["AR"] = subprocess.check_output(["which", "xiar"]).rstrip('\n')
    else:
        raise Exception("compiler_id not known")
    #this line only works is the build directory has been deleted
    cmake_compiler_args = shlex.split("-D CMAKE_C_COMPILER={} -D CMAKE_CXX_COMPILER={} "
                                      "-D CMAKE_LINKER={} -D CMAKE_AR={}"
                                      .format(env["CC"], env["CXX"], env["LD"], env["AR"]))
    return env, cmake_compiler_args

def run_cmake(compiler_id="unix"):
    if not os.path.isdir(cmake_build_dir):
        os.makedirs(cmake_build_dir)
    print "\nrunning cmake in directory", cmake_build_dir
    cwd = os.path.abspath(os.path.dirname(__file__))
    env, cmake_compiler_args = set_compiler_env(compiler_id)

    p = subprocess.call(["cmake"] + cmake_compiler_args + [cwd], cwd=cmake_build_dir, env=env)
    if p != 0:
        raise Exception("running cmake failed")
    print "\nbuilding files in cmake directory"
    if len(cmake_parallel_args) > 0:
        print "make flags:", cmake_parallel_args
    p = subprocess.call(["make"] + cmake_parallel_args, cwd=cmake_build_dir)
    if p != 0:
        raise Exception("building libraries with CMake Makefile failed")
    print "finished building the extension modules with cmake\n"


run_cmake(compiler_id=idcompiler)


# Now that the cython libraries are built, we have to make sure they are copied to
# the correct location.  This means in the source tree if build in-place, or
# somewhere in the build/ directory otherwise.  The standard distutils
# knows how to do this best.  We will overload the build_ext command class
# to simply copy the pre-compiled libraries into the right place
class build_ext_precompiled(old_build_ext):
    def build_extension(self, ext):
        """overload the function that build the extension

        This does nothing but copy the precompiled library stored in extension.sources[0]
        to the correct destination based on extension.name and whether it is an in-place build
        or not.
        """
        ext_path = self.get_ext_fullpath(ext.name)
        pre_compiled_library = ext.sources[0]
        if pre_compiled_library[-3:] != ".so":
            raise RuntimeError("library is not a .so file: " + pre_compiled_library)
        if not os.path.isfile(pre_compiled_library):
            raise RuntimeError("file does not exist: " + pre_compiled_library + " Did CMake not run correctly")
        print "copying", pre_compiled_library, "to", ext_path
        shutil.copy2(pre_compiled_library, ext_path)

# Construct extension modules for all the cxx files
# The `name` of the extension is, as usual, the python path (e.g. pele.optimize._lbfgs_cpp).
# The `source` of the extension is the location of the .so file
cxx_modules = []
for fname in cxx_files:
    name = fname.replace(".cxx", "")
    name = name.replace("/", ".")
    lname = os.path.basename(fname)
    lname = lname.replace(".cxx", ".so")
    pre_compiled_lib = os.path.join(cmake_build_dir, lname)
    cxx_modules.append(Extension(name, [pre_compiled_lib]))

setup(cmdclass=dict(build_ext=build_ext_precompiled),
      ext_modules=cxx_modules)
