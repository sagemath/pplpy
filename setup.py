#!/usr/bin/env python
from setuptools import setup
from distutils.command.build_ext import build_ext as _build_ext
from distutils.cmd import Command
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

with open('version.txt') as f:
    VERSION = f.read()[:-1]
    assert 2 <= len(VERSION.split('.')) <= 3

extensions = [
    Extension('ppl.linear_algebra',
        sources=['ppl/linear_algebra.pyx', 'ppl/ppl_shim.cc'],
        depends=['ppl/*']),
    Extension('ppl.mip_problem',
        sources=['ppl/mip_problem.pyx', 'ppl/ppl_shim.cc'],
        depends=['ppl/*']),
    Extension('ppl.polyhedron',
        sources = ['ppl/polyhedron.pyx', 'ppl/ppl_shim.cc'],
        depends=['ppl/*']),
    Extension('ppl.generator',
        sources = ['ppl/generator.pyx', 'ppl/ppl_shim.cc'],
        depends=['ppl/*']),
    Extension('ppl.constraint',
        sources = ['ppl/constraint.pyx', 'ppl/ppl_shim.cc'],
        depends=['ppl/*']),
    Extension('ppl.bit_arrays',
        sources = ['ppl/bit_arrays.pyx', 'ppl/ppl_shim.cc'],
        depends=['ppl/*'])
    ]

setup(
    name='pplpy',
    version=VERSION,
    description='Python PPL wrapper',
    long_description=open("README.rst").read(),
    author='Vincent Delecroix',
    author_email='vincent.delecroix@labri.fr',
    url='https://gitlab.com/videlec/pplpy',
    download_url='https://pypi.org/project/pplpy/#files',
    license='GPL v3',
    platforms=['any'],
    packages=['ppl'],
    package_dir={'ppl': 'ppl'},
    package_data={'ppl': ['*.pxd', '*.h', '*.hh']},
    install_requires=['Cython', 'cysignals', 'gmpy2'],  # For pip install, pip can't read setup_requires
    include_dirs=['ppl'] + sys.path,
    ext_modules=extensions,
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: C++",
        "Programming Language :: Python",
        "Development Status :: 3 - Alpha",
        "Operating System :: Unix",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    keywords=['polyhedron', 'polytope', 'convex', 'mathematics', 'ppl', 'milp', 'linear-programming'],
    cmdclass = {'build_ext': build_ext, 'test': TestCommand}
)
