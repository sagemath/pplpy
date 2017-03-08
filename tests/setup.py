from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import sys
import ppl

extensions = [
    Extension("testpplpy", ["testpplpy.pyx"],
        include_dirs=ppl.__path__,
        )
]

setup(
    name="test pplpy",
    ext_modules=cythonize(extensions, include_path=sys.path, language='c++')
)