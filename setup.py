from setuptools import setup, Extension, distutils, find_packages
from setuptools.command.build_ext import build_ext
import platform
import codecs
import sys
import os
import glob
import shutil

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

def find_cuda_config():
    cuda_home = os.environ.get("CUDA_HOME") or os.environ.get("CUDA_PATH")
    nvcc = None

    if cuda_home:
        candidate = os.path.join(cuda_home, "bin", "nvcc")
        if os.path.exists(candidate):
            nvcc = candidate

    if not nvcc:
        nvcc = shutil.which("nvcc")

    if not nvcc:
        return None

    if not cuda_home:
        cuda_home = os.path.dirname(os.path.dirname(nvcc))

    include_dir = os.path.join(cuda_home, "include")
    lib64_dir = os.path.join(cuda_home, "lib64")
    lib_dir = os.path.join(cuda_home, "lib")

    return {
        "home": cuda_home,
        "nvcc": nvcc,
        "include_dir": include_dir,
        "lib_dir": lib64_dir if os.path.exists(lib64_dir) else lib_dir,
    }


def setup_extension_modules(cuda_config=None):
    source_files = [
        "src/python.cpp",
        "src/molDataTable.cpp",
        "src/fpstore.cpp",
        "src/searchEngine.cpp",
        "src/utils.cpp",
        "src/h5_utils.cpp",
        "src/innerClustering.cpp",
    ]

    define_macros = []
    include_dirs = ["/src", get_pybind_include()]

    if cuda_config:
        source_files.append("src/innerClustering_cuda.cu")
        define_macros.append(("USE_CUDA", "1"))
        include_dirs.append(cuda_config["include_dir"])

    extension = Extension("uffpsim.uffpsimLib", sources=source_files, include_dirs=[ "/src", get_pybind_include() ], language="c++",)
    extension.include_dirs = include_dirs
    extension.define_macros = define_macros
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
        cuda_config = find_cuda_config()

        if cuda_config:
            if ".cu" not in self.compiler.src_extensions:
                self.compiler.src_extensions.append(".cu")

            original_compile = self.compiler._compile
            nvcc_path = cuda_config["nvcc"]

            def compile_with_nvcc(obj, src, ext, cc_args, extra_postargs, pp_opts):
                if src.endswith(".cu"):
                    postargs = []
                    if isinstance(extra_postargs, dict):
                        postargs = extra_postargs.get("nvcc", [])
                    self.compiler.spawn([nvcc_path, "-c", src, "-o", obj, "-Xcompiler", "-fPIC"] + pp_opts + postargs)
                else:
                    if isinstance(extra_postargs, dict):
                        extra_postargs = extra_postargs.get("cxx", [])
                    original_compile(obj, src, ext, cc_args, extra_postargs, pp_opts)

            self.compiler._compile = compile_with_nvcc

        f5_cflags, h5_ldflags = self.find_hdf5()
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])
        nvcc_opts = ["-O3", "--std=c++14", "-lineinfo"]

        if cuda_config and cuda_config["include_dir"]:
            opts.append("-I" + cuda_config["include_dir"])

        if cuda_config and cuda_config["lib_dir"]:
            link_opts.append("-L" + cuda_config["lib_dir"])
            link_opts.append("-lcudart")

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
            ext.define_macros.append(("VERSION_INFO", '"{}"'.format(self.distribution.get_version())))
            if cuda_config:
                ext.extra_compile_args = {
                    "cxx": opts + f5_cflags,
                    "nvcc": nvcc_opts + f5_cflags,
                }
            else:
                ext.extra_compile_args = opts + f5_cflags
            ext.extra_link_args = link_opts + h5_ldflags
        build_ext.build_extensions(self)


setup(
    packages=find_packages(),
    include_package_data=True,
    ext_modules=setup_extension_modules(find_cuda_config()),
    cmdclass={"build_ext": BuildExt},
)
