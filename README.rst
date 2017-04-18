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

TODO

- Needs explanation
- link to examples

Source
------

You can find the latest version of the source code on github:
https://github.com/videlec/pplpy/

Documentation
-------------

The documentation is available at http://pythonhosted.org/pplpy/

Requirements
------------

- `PPL <http://bugseng.com/products/ppl/>`_

- `gmp <https://gmplib.org/>`_

- `Cython <http://cython.org>`_

- `cygsignals <https://pypi.python.org/pypi/cysignals>`_

- `gmpy2 <https://pypi.python.org/pypi/gmpy2>`_

On Debian/Ubuntu systems these can be installed with::

    $ sudo apt-get install cython libgmp-dev libppl-dev
