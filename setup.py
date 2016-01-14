from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

setup(
    name='ppl',
    packages=['ppl'],
    ext_modules=[
        Extension('ppl',
            sources=['src/ppl.pyx', 'src/ppl_shim.cc'],
            libraries=['gmp','gmpxx','ppl','m'],
            language='c++'),

        Extension('pylong',
            sources=['src/pylong.pyx'],
            libraries=['gmp'],
            language='c'),
        ],
    cmdclass  = {'build_ext': build_ext}
)
