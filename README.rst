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

Source
------

You can find the latest version of the source code on github:
https://github.com/videlec/pplpy/


Requirements
------------

- PPL library with headers files

- `gmp <https://gmplib.org/>`_

- `Cython <http://cython.org>`_

On Debian/Ubuntu systems these can be installed with::

    $ sudo apt-get install cython libgmp-dev libppl-dev
