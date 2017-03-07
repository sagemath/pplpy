#!/usr/bin/env python
from setuptools import setup
from distutils.command.build_ext import build_ext as _build_ext
from setuptools.extension import Extension

import sys

# Adapted from Cython's new_build_ext
class build_ext(_build_ext):
    def finalize_options(self):
        # Check dependencies
        try:
            from Cython.Build.Dependencies import cythonize
        except ImportError as E:
            sys.stderr.write("Error: {0}\n".format(E))
            sys.stderr.write("The installation of ppl requires Cython\n")
            sys.exit(1)

        try:
            # We need the header files for cysignals at compile-time
            import cysignals
        except ImportError as E:
            sys.stderr.write("Error: {0}\n".format(E))
            sys.stderr.write("The installation of ppl requires cysignals\n")
            sys.exit(1)

        self.distribution.ext_modules[:] = cythonize(
            self.distribution.ext_modules, include_path=sys.path)
        _build_ext.finalize_options(self)

VERSION = open('version.txt').read()[:-1]

extensions = [
    Extension("ppl.cygmp.utils",
        sources = ["ppl/cygmp/utils.pyx"],
        depends = ["ppl/cygmp/*"]),
    Extension("ppl.cygmp.pylong",
        sources = ["ppl/cygmp/pylong.pyx"],
        depends = ["ppl/cygmp/*"]),
    Extension('ppl.wrappers',
        sources=['ppl/wrappers.pyx', 'ppl/ppl_shim.cc'],
        depends=['ppl/*'])
    ]

setup(
    name='pplpy',
    version=VERSION,
    description='Python PPL wrapper',
    long_description=open("README.rst").read(),
    author='Vincent Delecroix',
    author_email='vincent.delecroix@labri.fr',
    url='https://github.com/videlec/pplpy',
    download_url ='https://github.com/videlec/pplpy/archive/{}.tar.gz'.format(VERSION),
    license='GPL v3',
    platforms=['any'],
    packages=['ppl', 'ppl.cygmp'],
    package_dir={'ppl': 'ppl', 'ppl.cygmp': 'ppl/cygmp'},
    package_data={'ppl': ['*.pxd', '*.h', '*.hh'],
                  'ppl.cygmp': ['*.pxd', '*.h', '*.hh']},
    install_requires=['Cython', 'cysignals'],  # For pip install, pip can't read setup_requires
    include_dirs=['ppl'],
    ext_modules=extensions,
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: C++",
        "Programming Language :: Python",
        "Development Status :: 3 - Alpha",
        "Operating System :: Unix",
        "Intended Audience :: Science/Research",
    ],
    keywords=['polyhedron', 'polytope', 'convex', 'mathematics', 'ppl', 'milp', 'linear-programming'],
    cmdclass = {'build_ext': build_ext}
)
