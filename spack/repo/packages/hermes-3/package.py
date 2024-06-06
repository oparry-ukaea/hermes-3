from spack import *
import os
import shutil


class Hermes3(CMakePackage):
    """Hermes3"""

    homepage = "https://github.com/bendudson/hermes-3"

    git = "https://github.com/bendudson/hermes-3.git"
    version("working", branch="master")
    version("master", branch="master")

    depends_on("cmake@3.24:", type="build")
    #depends_on("adios2", type=("build", "link", "run"))
    depends_on("fftw", type=("build", "link", "run"))
    depends_on("mpi", type=("build", "link", "run"))
    depends_on("netcdf-cxx4", type=("build", "link", "run"))
    depends_on("py-cython", type=("build", "link", "run"))
    depends_on("py-jinja2", type=("build", "link", "run"))

    def cmake_args(self):
        args = []
        args.append("-DBOUT_DOWNLOAD_SUNDIALS=ON")
        return args
