PPL Python wrapper
==================

This Python package provides a wrapper to the C++ `Parma Polyhedra Library
(PPL) <http://bugseng.com/products/ppl/>`_.

The whole package is essentially a copy/paste of the corresponding file
from the `Sage <http://sagemath.org>`_ software. Though, the pplpy package
does not need Sage and can be used directly from Python.

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

Source
------

You can find the latest version of the source code on github:
https://github.com/videlec/pplpy/


Requirements
------------

- PPL libraries with headers

- `Cython <http://cython.org>`_
