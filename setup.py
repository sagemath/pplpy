#!/usr/bin/env python
from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext as _build_ext
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

        # Generate auto-generated sources from pari.desc
        #from autogen import rebuild
        #rebuild()

        self.distribution.ext_modules[:] = cythonize(
            self.distribution.ext_modules, include_path=sys.path)
        _build_ext.finalize_options(self)

VERSION = open('version.txt').read()[:-1]

extensions = [
    Extension("cygmp.utils",
        sources = ["src/cygmp/utils.pyx"],
        depends = ["src/cygmp/*"],
        libraries = ["gmp"],
        language = 'c'),
    Extension("cygmp.pylong",
        sources = ["src/cygmp/pylong.pyx"],
        depends = ["src/cygmp/*"],
        libraries = ["gmp"],
        language = 'c'),
    Extension('ppl',
        sources=['src/ppl.pyx', 'src/ppl_shim.cc'],
        depends=['src/ppl.pxd', 'src/ppl_shim.hh', 'src/ppl_decl.pxd'],
        libraries=['gmp','gmpxx','ppl','m'],
        language='c++'),
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
    packages=["ppl"],
    package_dir={'ppl': 'src'},
    package_data={'ppl': ['*.pxd', '*.h']},
    install_requires=['Cython', 'cysignals'],  # For pip install, pip can't read setup_requires
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
