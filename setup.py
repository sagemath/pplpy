#!/usr/bin/env python

import os
import sys

from setuptools import setup
from setuptools.config import read_configuration
from setuptools.extension import Extension

from distutils.command.build_ext import build_ext as _build_ext
from distutils.cmd import Command

dir_path = os.path.dirname(os.path.realpath(__file__))
conf_dict = read_configuration(os.path.join(dir_path, "setup.cfg"))
# NOTE: Python2.7 do not support multiple dictionaries unpacking
kwds = conf_dict['metadata']
kwds.update(conf_dict['options'])

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
            self.distribution.ext_modules,
            include_path=sys.path,
            compiler_directives={'embedsignature': True})

        _build_ext.finalize_options(self)

class TestCommand(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import subprocess, os, tempfile
        from shutil import rmtree
        from distutils.dir_util import copy_tree

        old_path = os.getcwd()
        tempdir_path = tempfile.mkdtemp()
        try:
            copy_tree('./tests', tempdir_path)
            os.chdir(tempdir_path)

            if subprocess.call([sys.executable, 'runtests.py']):
                raise SystemExit("Doctest failures")

            if subprocess.call([sys.executable, 'setup.py', 'build_ext', '--inplace']) or \
                    subprocess.call([sys.executable, '-c', "import testpplpy; testpplpy.test(); testpplpy.example()"]):
                raise SystemExit("Cython test 1 failure")

            if subprocess.call([sys.executable, 'setup2.py', 'build_ext', '--inplace']) or \
                    subprocess.call([sys.executable, '-c', "import testpplpy2; testpplpy2.test(); testpplpy2.example()"]):
                raise SystemExit("Cython test 2 failure")
        finally:
            os.chdir(old_path)
            rmtree(tempdir_path)

extensions = [
    Extension('ppl.linear_algebra', sources=['ppl/linear_algebra.pyx', 'ppl/ppl_shim.cc']),
    Extension('ppl.mip_problem', sources=['ppl/mip_problem.pyx', 'ppl/ppl_shim.cc']),
    Extension('ppl.polyhedron', sources = ['ppl/polyhedron.pyx', 'ppl/ppl_shim.cc']),
    Extension('ppl.generator', sources = ['ppl/generator.pyx', 'ppl/ppl_shim.cc']),
    Extension('ppl.constraint', sources = ['ppl/constraint.pyx', 'ppl/ppl_shim.cc']),
    Extension('ppl.congruence', sources=['ppl/congruence.pyx', 'ppl/ppl_shim.cc']),
    Extension('ppl.bit_arrays', sources = ['ppl/bit_arrays.pyx', 'ppl/ppl_shim.cc']),
    ]

setup(
    long_description = open("README.rst").read(),
    package_dir = {'ppl': 'ppl'},
    package_data = {'ppl': ['*.pxd', '*.h', '*.hh']},
    include_dirs = ['ppl'] + sys.path,
    ext_modules = extensions,
    cmdclass = {'build_ext': build_ext, 'test': TestCommand},
    **kwds
)
