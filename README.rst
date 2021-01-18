PPL Python wrapper
==================

This Python package provides a wrapper to the C++ `Parma Polyhedra Library
(PPL) <http://bugseng.com/products/ppl/>`_.

The whole package started as a fork of a tiny part of the `Sage
<http://sagemath.org>`_ software.

How it works
------------

The names of objects and methods are the same as in the library:

.. code:: python

    >>> import ppl
    >>> x = ppl.Variable(0)
    >>> y = ppl.Variable(1)
    >>> z = ppl.Variable(2)
    >>> cs = ppl.Constraint_System()
    >>> cs.insert(x >= 0)
    >>> cs.insert(y >= 0)
    >>> cs.insert(z >= 0)
    >>> cs.insert(x + y + z == 1)
    >>> poly = ppl.C_Polyhedron(cs)
    >>> poly.minimized_generators()
    Generator_System {point(1/1, 0/1, 0/1), point(0/1, 1/1, 0/1), point(0/1, 0/1, 1/1)}

The available objects and functions from the `ppl` Python module are:

- `Variable`, `Variables_Set`, `Linear_Expression` (defined in `ppl.linear_algebra`)

- `MIP_Problem` (defined in `ppl.mip_problem`)

- `C_Polyhedron`, `NNC_Polyhedron` (defined in `ppl.polyhedron`)

- `Generator`, `Generator_System`, `Poly_Gen_Relation`, `point`,
  `closure_point`, `ray`, `line` (defined in `ppl.generator`)

- `Constraint`, `Constraint_System`, `Poly_Con_Relation`,
  `inequality`, `equation`, `strict_inequality` (defined in `ppl.constraint`)

Installation
------------

The project is available at `Python Package Index <https://pypi.python.org/pypi/pplpy/>`_ and
can be installed with pip::

    $ pip install pplpy

Note that if you have gmp and ppl installed in a non standard directory (e.g. you use brew
on MacOSX) then you need to set appropriately the variables `CFLAGS` before calling `pip`. For
example::

    $ export CFLAGS="-I/path/to/gmp/include/ -L/path/to/gmp/lib/ -I/path/to/ppl/include/ -L/path/to/ppl/lib $CFLAGS"
    $ pip install pplpy

Using from Cython
-----------------

All Python classes from pplpy are extension types and can be used with Cython. Each
extension type carries an attribute `thisptr` that holds a pointer to
the corresponding C++ object from ppl.

A complete example is provided with the files `tests/testpplpy.pyx` and `tests/setup.py`.

Source
------

You can find the latest version of the source code on gitlab:
https://gitlab.com/videlec/pplpy

Documentation
-------------

An online version of the documentation is available at http://www.labri.fr/perso/vdelecro/pplpy/latest/

Compiling the html documentation requires make and `sphinx <http://www.sphinx-doc.org/en/master/>`_.
Before building the documentation, you need to install the pplpy package (sphinx uses Python introspection).
The documentation source code is contained in the repository `docs` where there is a standard
Makefile with a target `html`. Running `make html` in the `docs` repository builds the documentation
inside `docs/build/html`. For more configuration options, run `make help`.

License
-------

pplpy is distributed under the terms of the GNU General Public License (GPL)
published by the Free Software Foundation; either version 3 of
the License, or (at your option) any later version. See http://www.gnu.org/licenses/.

Requirements
------------

- `gmp <https://gmplib.org/>`_

- `mpfr <http://www.mpfr.org/>`_

- `mpc <http://www.multiprecision.org/index.php?prog=mpc>`_

- `PPL <http://bugseng.com/products/ppl/>`_

- `Cython <http://cython.org>`_

- `cysignals <https://pypi.python.org/pypi/cysignals>`_

- `gmpy2 <https://pypi.python.org/pypi/gmpy2>`_: version >= 2.1.0a4 from sources (see below)

On Debian/Ubuntu systems these can be installed with::

    $ sudo apt-get install libgmp-dev libmpfr-dev libmpc-dev libppl-dev cython
    $ pip install cysignals --user
    $ pip install gmpy2 --pre --user

The pip optional option `--user` allows to install python packages for a single
user with no need for administrator rights. The two pip install commands might
be replaced by `sudo pip install PKG` (not recommended). On recent Debian/Ubuntu systems,
cysignals is also available as a package under the name `python-cysignals` for
Python 2 and `python3-cysignals` for Python 3.
