from setuptools import setup
from Cython.Build import cythonize


opts = {
    'compiler_directives': {'language_level': '3'}
    }


setup(ext_modules=cythonize("testpplpy.pyx", **opts))
