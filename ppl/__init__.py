r"""
Cython wrapper for the Parma Polyhedra Library (PPL)

The Parma Polyhedra Library (PPL) is a library for polyhedral computations over
the rationals. This interface tries to reproduce the C++ API as faithfully as possible
in Python. For example, the following C++ excerpt:

.. code-block:: c++

    Variable x(0);
    Variable y(1);
    Constraint_System cs;
    cs.insert(x >= 0);
    cs.insert(x <= 3);
    cs.insert(y >= 0);
    cs.insert(y <= 3);
    C_Polyhedron poly_from_constraints(cs);

translates into:

>>> from ppl import Variable, Constraint_System, C_Polyhedron
>>> x = Variable(0)
>>> y = Variable(1)
>>> cs = Constraint_System()
>>> cs.insert(x >= 0)
>>> cs.insert(x <= 3)
>>> cs.insert(y >= 0)
>>> cs.insert(y <= 3)
>>> poly_from_constraints = C_Polyhedron(cs)

The same polyhedron constructed from generators:

>>> from ppl import Variable, Generator_System, C_Polyhedron, point
>>> gs = Generator_System()
>>> gs.insert(point(0*x + 0*y))
>>> gs.insert(point(0*x + 3*y))
>>> gs.insert(point(3*x + 0*y))
>>> gs.insert(point(3*x + 3*y))
>>> poly_from_generators = C_Polyhedron(gs)

Rich comparisons test equality/inequality and strict/non-strict
containment:

>>> poly_from_generators == poly_from_constraints
True
>>> poly_from_generators >= poly_from_constraints
True
>>> poly_from_generators <  poly_from_constraints
False
>>> poly_from_constraints.minimized_generators()
Generator_System {point(0/1, 0/1), point(0/1, 3/1), point(3/1, 0/1), point(3/1, 3/1)}
>>> poly_from_constraints.minimized_constraints()
Constraint_System {-x0+3>=0, -x1+3>=0, x0>=0, x1>=0}

As we see above, the library is generally easy to use. There are a few
pitfalls that are not entirely obvious without consulting the
documentation, in particular:

* There are no vectors used to describe :class:`Generator` (points,
  closure points, rays, lines) or :class:`Constraint` (strict
  inequalities, non-strict inequalities, or equations). Coordinates
  are always specified via linear polynomials in :class:`Variable`

* All coordinates of rays and lines as well as all coefficients of
  constraint relations are (arbitrary precision) integers. Only the
  generators :func:`point` and :func:`closure_point` allow one to
  specify an overall divisor of the otherwise integral
  coordinates. For example:

  >>> from ppl import Variable, point
  >>> x = Variable(0); y = Variable(1)
  >>> p = point( 2*x+3*y, 5 ); p
  point(2/5, 3/5)
  >>> p.coefficient(x)
  mpz(2)
  >>> p.coefficient(y)
  mpz(3)
  >>> p.divisor()
  mpz(5)

* PPL supports (topologically) closed polyhedra
  (:class:`C_Polyhedron`) as well as not necessarily closed polyhedra
  (:class:`NNC_Polyhedron`). Only the latter allows closure points
  (=points of the closure but not of the actual polyhedron) and strict
  inequalities (``>`` and ``<``)

The naming convention for the C++ classes is that they start with
``PPL_``, for example, the original ``Linear_Expression`` becomes
``PPL_Linear_Expression``. The Python wrapper has the same name as the
original library class, that is, just ``Linear_Expression``. In short:

* If you are using the Python wrapper (if in doubt: that's you), then
  you use the same names as the PPL C++ class library.

* If you are writing your own Cython code, you can access the
  underlying C++ classes by adding the prefix ``PPL_``.

Finally, PPL is fast. For example, here is the permutahedron of 5
basis vectors:

>>> from ppl import Variable, Generator_System, point, C_Polyhedron
>>> basis = range(0,5)
>>> x = [ Variable(i) for i in basis ]
>>> gs = Generator_System();
>>> from itertools import permutations
>>> for coeff in permutations(basis):
...    gs.insert(point( sum( (coeff[i]+1)*x[i] for i in basis ) ))
>>> C_Polyhedron(gs)
A 4-dimensional polyhedron in QQ^5 defined as the convex hull of 120 points

DIFFERENCES VS. C++

Since Python and C++ syntax are not always compatible, there are
necessarily some differences. The main ones are:

* The :class:`Linear_Expression` also accepts an iterable as input for
  the homogeneous coefficients.

AUTHORS:

- Volker Braun (2010): initial version (within Sage).
- Risan (2012): extension for MIP_Problem class (within Sage)
- Vincent Delecroix (2016-2020): convert Sage files into a standalone Python package,
  interface bit_array, interface congruence
- Vincent Klein (2017): improve doctest support and Python 3 compatibility
  Split the main code into several files.
  Remove the _mutable_immutable class.
"""

__version__ = "0.8.10rc2"

from .linear_algebra import (
        Variable, Variables_Set, Linear_Expression,
        )

from .mip_problem import MIP_Problem

from .polyhedron import C_Polyhedron, NNC_Polyhedron

from .generator import Generator, Generator_System, Poly_Gen_Relation
line = Generator.line
point = Generator.point
ray = Generator.ray
closure_point = Generator.closure_point

from .constraint import (
        Constraint, Constraint_System, Poly_Con_Relation,
        inequality, equation, strict_inequality)

from .congruence import Congruence, Congruence_System

from .bit_arrays import Bit_Row, Bit_Matrix
