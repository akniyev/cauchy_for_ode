from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize


examples_extension = Extension(
    name="laguerre_module",
    sources=["laguerre_module.pyx"],
    libraries=["laguerre"],
    library_dirs=["lib"],
    include_dirs=["lib"]
)

setup(
    name="laguerre_module",
    ext_modules=cythonize([examples_extension]), requires=['Cython']
)

