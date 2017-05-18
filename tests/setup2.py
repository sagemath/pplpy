from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import sys
import ppl

extensions = [
    Extension("testpplpy2", ["testpplpy2.pyx"],
        include_dirs=ppl.__path__ + sys.path,
        libraries=['ppl'],
        language='c++'
        )
]

setup(
    name="test pplpy2",
    ext_modules=cythonize(extensions, include_path=sys.path)
)
