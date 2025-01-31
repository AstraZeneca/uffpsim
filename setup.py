from setuptools import setup, Extension, distutils, find_packages
from setuptools.command.build_ext import build_ext
import platform
import codecs
import sys
import os
import glob

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked."""

    def __str__(self):
        import pybind11

        return pybind11.get_include()

def setup_extension_modules():
    source_files = ["src/python.cpp", "src/molDataTable.cpp", "src/fpstore.cpp", "src/searchEngine.cpp", "src/utils.cpp", "src/h5_utils.cpp", "src/innerClustering.cpp"]
    extension = Extension("uffpsim.uffpsimLib", sources=source_files, include_dirs=[ "/src", get_pybind_include() ], language="c++",)
    return [extension]

# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".cpp", delete=False) as f:
        f.write("int main (int argc, char **argv) { return 0; }")
        fname = f.name
    try:
        compiler.compile([fname], extra_postargs=[flagname])
    except distutils.errors.CompileError:
        return False
    finally:
        try:
            os.remove(fname)
        except OSError:
            pass
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14/17] compiler flag.
    The newer version is prefered over c++11 (when it is available).
    """
    flags = ["-std=c++17", "-std=c++14", "-std=c++11"]

    for flag in flags:
        if has_flag(compiler, flag):
            return flag

    raise RuntimeError("Unsupported compiler -- at least C++11 support " "is needed!")


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""

    c_opts = {"msvc": ["/EHsc", "/arch:AVX"], "unix": ["-O3", "-funroll-all-loops", "-fopenmp"]}

    l_opts = {"msvc": [], "unix": ["-fopenmp"]}

    if sys.platform == "darwin":
        darwin_opts = ["-stdlib=libc++", "-mmacosx-version-min=10.9"]
        c_opts["unix"] += darwin_opts
        l_opts["unix"] += darwin_opts

    def find_hdf5(self):
        if 'CONDA_PREFIX' in os.environ:
            conda_path = os.environ["CONDA_PREFIX"]
            if os.path.exists(os.path.join(conda_path, "include", "H5Cpp.h")) and os.path.exists(os.path.join(conda_path, "lib", "libhdf5_cpp.so")):
                return ["-I" + os.path.join(conda_path, "include")], ["-L" + os.path.join(conda_path, "lib"), "-lhdf5_cpp", "-lhdf5_hl"]

        try:
            import pkgconfig
        except ImportError:
            pkgconfig = None
            
        if not pkgconfig:
            raise LookupError("pkgconfig not found! cannot search for HDF5 dev library. Install pkgconfig and try again.")
        
        if not pkgconfig.exists("hdf5"):
            raise LookupError("HDF5 dev library not found")

        cflags = pkgconfig.cflags("hdf5")
        ldflags = pkgconfig.libs("hdf5")

        if not cflags or not ldflags:
            raise LookupError("Could not find HDF5 dev library through pkgconfig!")

        return cflags.split(), ldflags.split() + ['-lhdf5_cpp', '-lhdf5_hl']

    def build_extensions(self):
        f5_cflags, h5_ldflags = self.find_hdf5()
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])
        if ct == "unix":
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, "-fvisibility=hidden"):
                opts.append("-fvisibility=hidden")

            if 'UFFPSIM_NATIVE' in os.environ:
                opts.append("-march=native")

            elif 'UFFPSIM_AVX512' in os.environ and has_flag(self.compiler, "-mavx512vpopcntdq"):
                opts.append("-DHAVE_AVX512VPOPCNTDQ")
                opts.append("-mavx512vpopcntdq")

            else:
                opts.append("-msse4.2")

        for ext in self.extensions:
            ext.define_macros = [
                ("VERSION_INFO", '"{}"'.format(self.distribution.get_version()))
            ]
            ext.extra_compile_args = opts + f5_cflags
            ext.extra_link_args = link_opts + h5_ldflags
        build_ext.build_extensions(self)


setup(
    packages=find_packages(),
    ext_modules=setup_extension_modules(),
    cmdclass={"build_ext": BuildExt},
)
