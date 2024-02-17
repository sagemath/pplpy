.. pplpy documentation master file

Welcome to pplpy's documentation!
=================================

Installation
------------

pplpy is available from the `Python Package Index <https://pypi.org/project/pplpy/>`_. You can hence install it with::

    $ pip install pplpy

The code source can be obtained from `github <https://github.com/sagemath/pplpy>`_::

    $ git clone https://github.com/sagemath/pplpy.git

Introduction
------------

.. automodule:: ppl

Using in Cython
---------------

All types from ppl are extension types and can be used with Cython. The following is a
short sample of Cython code using pplpy::

    from ppl.linear_algebra cimport Variable
    from ppl.constraint cimport Constraint_System
    from ppl.polyhedron cimport C_Polyhedron

    cdef Variable x = ppl.Variable(0)
    cdef Variable y = ppl.Variable(1)
    cdef Variable z = ppl.Variable(2)
    cdef Constraint_System cs = Constraint_System()
    cs.insert(x >= 0)
    cs.insert(y >= 0)
    cs.insert(x + y + z == 1)
    cdef C_Polyhedron poly = C_Polyhedron(cs)
    print(poly.minimized_generators())

Each extension type carries an attribute ``thisptr`` that holds a pointer to
the corresponding C++ object from ppl. Continuing the above example, one can
do::

    print('dim = %lu' % poly.thisptr.space_dimension())

To avoid name collisions, the original C++ class names are prefixed with
``PPL_``, for example, the original ``Linear_Expression`` becomes
``PPL_Linear_Expression``. The Python wrapper has the same name as the original
library class, that is, just ``Linear_Expression``. All ``ppl`` declarations
are done in the ``.pxd`` file ``ppl/ppl_decl.pxd``. Only few functionalities
from ``ppl`` are exposed in ``ppl_decl.pxd``. It is also preferable to avoid
mixing C++ ``ppl`` objects with Cython ``pplpy`` extension types.

To compile a Cython extension using pplpy you need to:

- add the path of your pplpy installation to ``include_dirs``
- link with the ``ppl`` library

Here is a minimal example of a ``setup.py`` file for a unique Cython file
called ``test.pyx``::

    from distutils.core import setup
    from distutils.extension import Extension
    from Cython.Build import cythonize
    import ppl

    extensions = [
        Extension("test", ["test.pyx"], libraries=['ppl'], include_dirs=ppl.__path__, language='c++')
    ]

    setup(ext_modules=cythonize(extensions))

Module `ppl.constraint`
-----------------------

.. automodule:: ppl.constraint
    :members:
    :undoc-members:
    :show-inheritance:

Module `ppl.linear_algebra`
---------------------------

.. automodule:: ppl.linear_algebra
    :members:
    :undoc-members:
    :show-inheritance:

Module `ppl.generator`
----------------------

.. automodule:: ppl.generator
    :members:
    :undoc-members:
    :show-inheritance:

Module `ppl.polyhedron`
-----------------------

.. automodule:: ppl.polyhedron
    :members:
    :undoc-members:
    :show-inheritance:

Module `ppl.mip_problem`
------------------------

.. automodule:: ppl.mip_problem
    :members:
    :undoc-members:
    :show-inheritance:

Module `ppl.congruence`
-----------------------

.. automodule:: ppl.congruence
    :members:
    :undoc-members:
    :show-inheritance:
