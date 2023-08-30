from setuptools import Extension, setup
from Cython.Build import cythonize


extension = Extension("testpplpy2", ["testpplpy2.pyx"], language='c++', libraries=['ppl'])
opts = {
    'compiler_directives': {'language_level': '3'}
    }


setup(ext_modules=cythonize(extension, *opts))
