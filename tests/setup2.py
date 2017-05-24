from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import ppl

extensions = [
    Extension("testpplpy2", ["testpplpy2.pyx"],
        include_dirs=ppl.__path__,
        libraries=['ppl'],
        language='c++'
        )
]

setup(ext_modules=cythonize(extensions))
