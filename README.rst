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

You can find the latest version of the source code on github:
https://github.com/videlec/pplpy/

Documentation
-------------

The documentation is available at http://pythonhosted.org/pplpy/

Requirements
------------

- `gmp <https://gmplib.org/>`_

- `mpfr <http://www.mpfr.org/>`_

- `mpc <http://www.multiprecision.org/index.php?prog=mpc>`_

- `PPL <http://bugseng.com/products/ppl/>`_

- `Cython <http://cython.org>`_

- `cygsignals <https://pypi.python.org/pypi/cysignals>`_

- `gmpy2 <https://pypi.python.org/pypi/gmpy2>`_

On Debian/Ubuntu systems these can be installed with::

    $ sudo apt-get install libgmp-dev libmpfr-dev libmpc-dev libppl-dev cython
    $ sudo pip install cysignals [--user]
    $ pip install gmpy2==2.1.0a1 --no-binary ":all:" [--user]

The pip optional option `--user` allows to install python packages for a single
user with no need for administrator rights.
