import os.path
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import ppl


extension = Extension("testpplpy", ["testpplpy.pyx"])

opts = {
    'compiler_directives': {'language_level': '3'}
    }

setup(ext_modules=cythonize(extension, **opts))
