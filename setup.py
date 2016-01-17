#!/usr/bin/env python
from setuptools import setup
from setuptools.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

VERSION = open('version.txt').read()[:-1]

setup(
    name='pplpy',
    version=VERSION,
    description='Python PPL wrapper',
    long_description=open("README.rst").read(),
    author='Vincent Delecroix',
    author_email='vincent.delecroix@labri.fr',
    url = 'https://github.com/videlec/pplpy',
    download_url ='https://github.com/videlec/pplpy/archive/{}.tar.gz'.format(VERSION),
    license='GPL v3',
    platforms=['any'],
    ext_modules=[
        Extension('ppl',
            sources=['src/ppl.pyx', 'src/ppl_shim.cc'],
            depends=['src/ppl.pxd', 'src/ppl_shim.hh', 'src/ppl_decl.pxd'],
            libraries=['gmp','gmpxx','ppl','m'],
            language='c++'),

        Extension('pylong',
            sources=['src/pylong.pyx'],
            depends=['src/pylong.pxd'],
            libraries=['gmp'],
            language='c'),
        ],
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: C++",
        "Programming Language :: Python",
        "Development Status :: 3 - Alpha",
        "Operating System :: Unix",
        "Intended Audience :: Science/Research",
    ],
    keywords=['polyhedron', 'polytope', 'convex', 'mathematics', 'ppl', 'milp', 'linear-programming'],
    cmdclass  = {'build_ext': build_ext}
)
