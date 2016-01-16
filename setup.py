#!/usr/bin/env python
from setuptools import setup
from setuptools.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

setup(
    name='pplpy',
    version='0.1',
    description='Python PPL wrapper',
    long_description=open("README.txt").read(),
    author='Vincent Delecroix',
    author_email='vincent.delecroix@labri.fr',
    url = 'https://github.com/videlec/pplpy',
    download_url = 'https://github.com/videlec/pplpy/tarball/0.1',
    license='GPL v3',
    platforms=['any'],
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
    classifiers=[
        "License :: OSI Approved :: GNU General Public License (GPL) version 3",
        "Programming Language :: Python",
        "Development Status :: 0 - Beta",
        "Operating System :: Unix",
    ],
    keywords=['polyhedron', 'polytope', 'convex', 'mathematics', 'ppl', 'milp', 'linear-programming'],
    cmdclass  = {'build_ext': build_ext}
)
