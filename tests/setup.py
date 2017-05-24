from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import ppl

extensions = [
    Extension("testpplpy", ["testpplpy.pyx"], include_dirs=ppl.__path__)
]

setup(ext_modules=cythonize(extensions))
