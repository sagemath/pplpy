from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

setup(
    name='ppl',
    version='0.0',
    description='PPL wrapper',
    author='Vincent Delecroix',
    author_email='vincent.delecroix@labri.fr',
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
