# distutils: language = c++
# distutils: libraries = gmp gmpxx ppl m
r"""
Cython wrapper for the Parma Polyhedra Library (PPL)

The Parma Polyhedra Library (PPL) is a library for polyhedral computations over
`\QQ`. This interface tries to reproduce the C++ API as faithfully as possible
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
  2
  >>> p.coefficient(y)
  3
  >>> p.divisor()
  5

* PPL supports (topologically) closed polyhedra
  (:class:`C_Polyhedron`) as well as not neccesarily closed polyhedra
  (:class:`NNC_Polyhedron`). Only the latter allows closure points
  (=points of the closure but not of the actual polyhedron) and strict
  inequalities (``>`` and ``<``)

The naming convention for the C++ classes is that they start with
``PPL_``, for example, the original ``Linear_Expression`` becomes
``PPL_Linear_Expression``. The Python wrapper has the same name as the
original library class, that is, just ``Linear_Expression``. In short:

* If you are using the Python wrapper (if in doubt: thats you), then
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
  the homogeneous cooefficients.

* :class:`Polyhedron` and its subclasses as well as
  :class:`Generator_System` and :class:`Constraint_System` can be set
  immutable via a ``set_immutable()`` method. This is the analog of
  declaring a C++ instance ``const``. All other classes are immutable
  by themselves.

AUTHORS:

- Volker Braun (2010-10-08): initial version (within Sage).
- Risan (2012-02-19): extension for MIP_Problem class (within Sage)
- Vincent Delecroix (2016): convert Sage files into a standalone Python package
"""

#*****************************************************************************
#       Copyright (C) 2010 Volker Braun  <vbraun.name@gmail.com>
#                     2016 Vincent Delecroix <vincent.delecroix@labri.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from cpython.int cimport PyInt_CheckExact
from cpython.long cimport PyLong_CheckExact

from pylong cimport mpz_get_pyintlong, mpz_set_pylong

try:
    from sage.all import Rational
    def Fraction(p,q):
        return Rational((p,q))
except ImportError:
    from fractions import Fraction

#TODO: find out whether it would be easy to switch between mpir and gmp
# these are just gmp declarations
# TODO:
# here we need to allow to branch this code with integer and rational classes.
# Give three options:
#  1) with Integer/Rational from Sage as below
#  2) with gmpy2
#  3) with Python native types
# and potentially more if the user want to dig in
#from sage.rings.integer cimport Integer
#from sage.rings.rational cimport Rational

# TODO: interruption buisness. This is internal to Sage. Though by default
# we could map sig_on/sig_off to a no-op
# how can these be changed at Python launched time? -> function pointers!
#include 'sage/ext/interrupt.pxi'


####################################################
# Potentially expensive operations:
#  - compute dual description
#  - solve linear program
# These can only be triggered by methods in the Polyhedron class
# they need to be wrapped in sig_on() / sig_off()
####################################################

cdef sig_on():
    pass
cdef sig_off():
    pass

####################################################
# PPL can use floating-point arithmetic to compute integers
cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library":
    cdef void set_rounding_for_PPL()
    cdef void restore_pre_PPL_rounding()

# but with PPL's rounding the gsl will be very unhappy; must turn off!
restore_pre_PPL_rounding()

cdef PPL_Coefficient PPL_Coefficient_from_pyobject(c):
    cdef mpz_t coeff
    if PyInt_CheckExact(c):
        return PPL_Coefficient(<long> c)
    elif PyLong_CheckExact(c):
        mpz_init(coeff)
        mpz_set_pylong(coeff, c)
        return PPL_Coefficient(coeff)
    else:
        try:
            return PPL_Coefficient_from_pyobject(c.__long__())
        except AttributeError:
            raise ValueError("unknown input type {}".format(type(c)))

cdef PPL_GeneratorType_str(PPL_GeneratorType t):
    if t == LINE:
        return 'line'
    elif t == RAY:
        return 'ray'
    elif t == POINT:
        return 'point'
    elif t == CLOSURE_POINT:
        return 'closure_point'
    else:
        raise RuntimeError

cdef _wrap_Constraint_System(PPL_Constraint_System constraint_system):
    """
    Wrap a C++ ``PPL_Constraint_System`` into a Cython ``Constraint_System``.
    """
    cdef Constraint_System cs = Constraint_System()
    del cs.thisptr
    cs.thisptr = new PPL_Constraint_System(constraint_system)
    return cs

cdef _wrap_Poly_Gen_Relation(PPL_Poly_Gen_Relation relation):
    """
    Wrap a C++ ``PPL_Poly_Gen_Relation`` into a Cython ``Poly_Gen_Relation``.
    """
    cdef Poly_Gen_Relation rel = Poly_Gen_Relation(True)
    rel.thisptr = new PPL_Poly_Gen_Relation(relation)
    return rel

cdef _wrap_Poly_Con_Relation(PPL_Poly_Con_Relation relation):
    """
    Wrap a C++ ``PPL_Poly_Con_Relation`` into a Cython ``Poly_Con_Relation``.
    """
    cdef Poly_Con_Relation rel = Poly_Con_Relation(True)
    rel.thisptr = new PPL_Poly_Con_Relation(relation)
    return rel

cdef _wrap_Generator_System(PPL_Generator_System generator_system):
    """
    Wrap a C++ ``PPL_Generator_System`` into a Cython ``Generator_System``.
    """
    cdef Generator_System gs = Generator_System()
    del gs.thisptr
    gs.thisptr = new PPL_Generator_System(generator_system)
    return gs

cdef _wrap_Constraint(PPL_Constraint constraint):
    """
    Wrap a C++ ``PPL_Constraint`` into a Cython ``Constraint``.
    """
    cdef Constraint c = Constraint(True)
    c.thisptr = new PPL_Constraint(constraint)
    return c

####################################################
### _mutable_or_immutable ##########################
####################################################
# TODO: remove it?
cdef class _mutable_or_immutable(object):
    r"""
    A base class for mutable or immutable objects.

    By default, any object is mutable. It can then be
    :meth:`set_immutable`.

    Examples:

        >>> from ppl import _mutable_or_immutable as ExampleObj
        >>> x = ExampleObj()
        >>> x.is_mutable()
        True
        >>> x.is_immutable()
        False
        >>> x.set_immutable()
        >>> x.is_mutable()
        False
    """
    def __cinit__(self):
        """
        The Cython constructor.

        Tests:

            >>> from ppl import _mutable_or_immutable as ExampleObj
            >>> x = ExampleObj()    # indirect doctest
            >>> x.is_mutable()
            True
        """
        self._is_mutable = True

    def set_immutable(self):
        """
        Make this object immutable.

        This operation cannot be undone.

        Examples:

        >>> from ppl import _mutable_or_immutable as ExampleObj
        >>> x = ExampleObj()
        >>> x.is_mutable()
        True
        >>> x.set_immutable()
        >>> x.is_mutable()
        False
        """
        self._is_mutable = False

    def is_mutable(self):
        """
        Return whether this object is mutable.

        The data members of the object can only be modified if the
        object is mutable.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import _mutable_or_immutable as ExampleObj
        >>> x = ExampleObj()
        >>> x.is_mutable()
        True
        """
        return self._is_mutable

    def is_immutable(self):
        """
        Return whether this object is immutable.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import _mutable_or_immutable as ExampleObj
        >>> x = ExampleObj()
        >>> x.is_immutable()
        False
        """
        return not self._is_mutable

    def assert_mutable(self, msg):
        r"""
        Raise ``ValueError`` if the object is not mutable.

        INPUT:

        - ``msg`` -- a string. The message to be returned together
          with the ``ValueError``

        OUTPUT:

        This method returns no output. A ``ValueError``` is raised if
        the object is not mutable.

        Examples:

        >>> from ppl import _mutable_or_immutable as ExampleObj
        >>> x = ExampleObj()
        >>> x.assert_mutable("this will not trigger")
        >>> x.set_immutable()
        >>> x.assert_mutable("this will trigger")
        Traceback (most recent call last):
        ...
        ValueError: this will trigger
        """
        if not self._is_mutable:
            raise ValueError(msg)

####################################################
### MIP_Problem ####################################
####################################################
cdef class MIP_Problem(object):
    r"""
    wrapper for PPL's MIP_Problem class

    An object of the class MIP_Problem represents a Mixed Integer
    (Linear) Program problem.

    INPUT:

    - ``dim`` -- integer
    - ``args`` -- an array of the defining data of the MIP_Problem.
      For each element, any one of the following is accepted:

      * A :class:`Constraint_System`.

      * A :class:`Linear_Expression`.

    OUTPUT:

    A :class:`MIP_Problem`.

    Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert(x >= 0)
        >>> cs.insert(y >= 0)
        >>> cs.insert(3 * x + 5 * y <= 10)
        >>> m = MIP_Problem(2, cs, x + y)
        >>> m
        MIP Problem (maximization, 2 variables, 3 constraints)
        >>> m.optimal_value()
        Fraction(10, 3)
        >>> float(_)
        3.3333333333333335
        >>> m.optimizing_point()
        point(10/3, 0/3)
    """
    def __repr__(self):
        """
        String representation of MIP Problem.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0 )
        >>> cs.insert( y >= 0 )
        >>> cs.insert( 3 * x + 5 * y <= 10 )
        >>> m = MIP_Problem(2, cs, x + y)
        >>> m
        MIP Problem (maximization, 2 variables, 3 constraints)
        """
        dim = self.space_dimension()
        ncs = sum(1 for _ in self)
        return 'MIP Problem ({}, {} variable{}, {} constraint{})'.format(
                self.optimization_mode(),
                dim,
                's' if dim > 1 else '',
                ncs,
                's' if ncs > 1 else '')

    def __cinit__(self, PPL_dimension_type dim = 0, *args):
        """
        The Cython constructor.

        Tests:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> MIP_Problem(0)
        A MIP_Problem
        Maximize: 0
        Subject to constraints
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert(x + y <= 2)
        >>> M = MIP_Problem(2, cs, 0)
        >>> M = MIP_Problem(2, cs, x)
        >>> M = MIP_Problem(2, None, None)
        Traceback (most recent call last):
        ...
        TypeError: Cannot convert NoneType to sage.libs.ppl.Constraint_System
        >>> M = MIP_Problem(2, cs, 'hey')
        Traceback (most recent call last):
        ...
        TypeError: unable to convert 'hey' to an integer
        >>> M = MIP_Problem(2, cs, x, 'middle')
        Traceback (most recent call last):
        ...
        ValueError: unknown mode 'middle'
        """
        if not args:
            self.thisptr = new PPL_MIP_Problem(dim)
            return
        elif len(args) == 1:
            raise ValueError('cannot initialize from {}'.format(args))

        cdef Constraint_System cs = <Constraint_System?>args[0]
        cdef Linear_Expression obj
        try:
            obj = <Linear_Expression?> args[1]
        except TypeError:
            obj = Linear_Expression(args[1])

        cdef PPL_Optimization_Mode mode = MAXIMIZATION
        if len(args) == 3:
            if args[2] == 'maximization':
                mode = MAXIMIZATION
            elif args[2] == 'minimization':
                mode = MINIMIZATION
            else:
                raise ValueError('unknown mode {!r}'.format(args[2]))

        self.thisptr = new PPL_MIP_Problem(dim, cs.thisptr[0], obj.thisptr[0], mode)

    def __dealloc__(self):
        """
        The Cython destructor
        """
        del self.thisptr

    def __iter__(self):
        r"""
        Iterator through the constraints

        Tests:

        >>> from ppl import Variable, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> M = MIP_Problem(2)
        >>> for c in M: print c
        >>> M.add_constraint(x + y <= 5)
        >>> for c in M: print c
        -x0-x1+5>=0
        >>> M.add_constraint(3*x - 18*y >= -2)
        >>> for c in M: print c
        -x0-x1+5>=0
        3*x0-18*x1+2>=0
        """
        return MIP_Problem_constraints_iterator(self)

    def constraints(self):
        r"""
        Return the constraints of this MIP

        The output is an instance of :class:`Constraint_System`.

        Examples:

        >>> from ppl import Variable, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> M = MIP_Problem(2)
        >>> M.add_constraint(x + y <= 5)
        >>> M.add_constraint(3*x - 18*y >= -2)
        >>> M.constraints()
        Constraint_System {-x0-x1+5>=0, 3*x0-18*x1+2>=0}

        Note that modifying the output of this method will not modify the
        underlying MIP problem object:

        >>> cs = M.constraints()
        >>> cs.insert(x <= 3)
        >>> cs
        Constraint_System {-x0-x1+5>=0, 3*x0-18*x1+2>=0, -x0+3>=0}
        >>> M.constraints()
        Constraint_System {-x0-x1+5>=0, 3*x0-18*x1+2>=0}
        """
        cdef Constraint_System c = Constraint_System(None)
        c.thisptr = mip_constraints(self.thisptr[0])
        return c

    def optimization_mode(self):
        """
        Return the optimization mode used in the MIP_Problem.

        It will return "maximization" if the MIP_Problem was set
        to MAXIMIZATION mode, and "minimization" otherwise.

        Examples:

            >>> from ppl import MIP_Problem
            >>> m = MIP_Problem()
            >>> m.optimization_mode()
            'maximization'
        """
        if self.thisptr.optimization_mode() == MAXIMIZATION:
            return "maximization"
        elif self.thisptr.optimization_mode() == MINIMIZATION:
            return "minimization"

    def optimal_value(self):
        """
        Return the optimal value of the MIP_Problem. ValueError thrown if self does not
        have an optimizing point, i.e., if the MIP problem is unbounded or not satisfiable.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0 )
        >>> cs.insert( y >= 0 )
        >>> cs.insert( 3 * x + 5 * y <= 10 )
        >>> m = MIP_Problem(2, cs, x + y)
        >>> m.optimal_value()
        Fraction(10, 3)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0 )
        >>> m = MIP_Problem(1, cs, x + x )
        >>> m.optimal_value()
        Traceback (most recent call last):
        ...
        ValueError: PPL::MIP_Problem::optimizing_point():
        *this does not have an optimizing point.
        """
        cdef PPL_Coefficient sup_n
        cdef PPL_Coefficient sup_d

        sig_on()
        try:
            self.thisptr.optimal_value(sup_n, sup_d)
        finally:
            sig_off()

        return Fraction(mpz_get_pyintlong(sup_n.get_mpz_t()), mpz_get_pyintlong(sup_d.get_mpz_t()))

    def space_dimension(self):
        """
        Return the space dimension of the MIP_Problem.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0)
        >>> cs.insert( y >= 0 )
        >>> cs.insert( 3 * x + 5 * y <= 10 )
        >>> m = MIP_Problem(2, cs, x + y)
        >>> m.space_dimension()
        2
        """
        return self.thisptr.space_dimension()

    def objective_function(self):
        """
        Return the optimal value of the MIP_Problem.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0)
        >>> cs.insert( y >= 0 )
        >>> cs.insert( 3 * x + 5 * y <= 10 )
        >>> m = MIP_Problem(2, cs, x + y)
        >>> m.objective_function()
        x0+x1
        """
        rc = Linear_Expression()
        rc.thisptr[0] = self.thisptr.objective_function()
        return rc

    def clear(self):
        """
        Reset the MIP_Problem to be equal to the trivial MIP_Problem.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0)
        >>> cs.insert( y >= 0 )
        >>> cs.insert( 3 * x + 5 * y <= 10 )
        >>> m = MIP_Problem(2, cs, x + y)
        >>> m.objective_function()
        x0+x1
        >>> m.clear()
        >>> m.objective_function()
        0
        """
        self.thisptr.clear()

    def add_space_dimensions_and_embed(self, PPL_dimension_type m):
        """
        Adds m new space dimensions and embeds the old MIP problem in the new vector space.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0)
        >>> cs.insert( y >= 0 )
        >>> cs.insert( 3 * x + 5 * y <= 10 )
        >>> m = MIP_Problem(2, cs, x + y)
        >>> m.add_space_dimensions_and_embed(5)
        >>> m.space_dimension()
        7
        """
        sig_on()
        self.thisptr.add_space_dimensions_and_embed(m)
        sig_off()

    def add_constraint(self, Constraint c):
        """
        Adds a copy of constraint c to the MIP problem.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> m = MIP_Problem()
        >>> m.add_space_dimensions_and_embed(2)
        >>> m.add_constraint(x >= 0)
        >>> m.add_constraint(y >= 0)
        >>> m.add_constraint(3 * x + 5 * y <= 10)
        >>> m.set_objective_function(x + y)
        >>> m.optimal_value()
        Fraction(10, 3)

        Tests:

        >>> z = Variable(2)
        >>> m.add_constraint(z >= -3)
        Traceback (most recent call last):
        ...
        ValueError: PPL::MIP_Problem::add_constraint(c):
        c.space_dimension() == 3 exceeds this->space_dimension == 2.
        """
        sig_on()
        try:
            self.thisptr.add_constraint(c.thisptr[0])
        finally:
            sig_off()

    def add_constraints(self, Constraint_System cs):
        """
        Adds a copy of the constraints in cs to the MIP problem.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert(x >= 0)
        >>> cs.insert(y >= 0)
        >>> cs.insert(3 * x + 5 * y <= 10)
        >>> m = MIP_Problem(2)
        >>> m.set_objective_function(x + y)
        >>> m.add_constraints(cs)
        >>> m.optimal_value()
        Fraction(10, 3)

        Tests:

        >>> p = Variable(9)
        >>> cs.insert(p >= -3)
        >>> m.add_constraints(cs)
        Traceback (most recent call last):
        ...
        ValueError: PPL::MIP_Problem::add_constraints(cs):
        cs.space_dimension() == 10 exceeds this->space_dimension() == 2.
        """
        sig_on()
        try:
            self.thisptr.add_constraints(cs.thisptr[0])
        finally:
            sig_off()

    def set_objective_function(self, Linear_Expression obj):
        """
        Sets the objective function to obj.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> m = MIP_Problem()
        >>> m.add_space_dimensions_and_embed(2)
        >>> m.add_constraint(x >= 0)
        >>> m.add_constraint(y >= 0)
        >>> m.add_constraint(3 * x + 5 * y <= 10)
        >>> m.set_objective_function(x + y)
        >>> m.optimal_value()
        Fraction(10, 3)

        Tests:

        >>> z = Variable(2)
        >>> m.set_objective_function(x + y + z)
        Traceback (most recent call last):
        ...
        ValueError: PPL::MIP_Problem::set_objective_function(obj):
        obj.space_dimension() == 3 exceeds this->space_dimension == 2.
        """
        self.thisptr.set_objective_function(obj.thisptr[0])

    def set_optimization_mode(self, mode):
        """
        Sets the optimization mode to mode.

        Examples:

        >>> from ppl import MIP_Problem
        >>> m = MIP_Problem()
        >>> m.optimization_mode()
        'maximization'
        >>> m.set_optimization_mode('minimization')
        >>> m.optimization_mode()
        'minimization'

        Tests:

        >>> m.set_optimization_mode('max')
        Traceback (most recent call last):
        ...
        ValueError: Unknown value: mode=max.
        """
        if mode == 'minimization':
            self.thisptr.set_optimization_mode(MINIMIZATION)
        elif mode == 'maximization':
            self.thisptr.set_optimization_mode(MAXIMIZATION)
        else:
            raise ValueError, 'Unknown value: mode='+str(mode)+'.'

    def is_satisfiable(self):
        """
        Check if the MIP_Problem is satisfiable

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> m = MIP_Problem()
        >>> m.add_space_dimensions_and_embed(2)
        >>> m.add_constraint(x >= 0)
        >>> m.add_constraint(y >= 0)
        >>> m.add_constraint(3 * x + 5 * y <= 10)
        >>> m.is_satisfiable()
        True
        """
        return self.thisptr.is_satisfiable()

    def evaluate_objective_function(self, Generator evaluating_point):
        """
        Return the result of evaluating the objective function on evaluating_point. ValueError thrown
        if self and evaluating_point are dimension-incompatible or if the generator evaluating_point is not a point.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem, Generator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> m = MIP_Problem()
        >>> m.add_space_dimensions_and_embed(2)
        >>> m.add_constraint(x >= 0)
        >>> m.add_constraint(y >= 0)
        >>> m.add_constraint(3 * x + 5 * y <= 10)
        >>> m.set_objective_function(x + y)
        >>> g = Generator.point(5 * x - 2 * y, 7)
        >>> m.evaluate_objective_function(g)
        Fraction(3, 7)
        >>> z = Variable(2)
        >>> g = Generator.point(5 * x - 2 * z, 7)
        >>> m.evaluate_objective_function(g)
        Traceback (most recent call last):
        ...
        ValueError: PPL::MIP_Problem::evaluate_objective_function(p, n, d):
        *this and p are dimension incompatible.
        """
        cdef PPL_Coefficient sup_n
        cdef PPL_Coefficient sup_d

        sig_on()
        try:
            self.thisptr.evaluate_objective_function(evaluating_point.thisptr[0], sup_n, sup_d)
        finally:
            sig_off()

        return Fraction(mpz_get_pyintlong(sup_n.get_mpz_t()),
                mpz_get_pyintlong(sup_d.get_mpz_t()))

    def solve(self):
        """
        Optimizes the MIP_Problem

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> m = MIP_Problem()
        >>> m.add_space_dimensions_and_embed(2)
        >>> m.add_constraint(x >= 0)
        >>> m.add_constraint(y >= 0)
        >>> m.add_constraint(3 * x + 5 * y <= 10)
        >>> m.set_objective_function(x + y)
        >>> m.solve()
        {'status': 'optimized'}
        """
        sig_on()
        try:
            tmp = self.thisptr.solve()
        finally:
            sig_off()
        if tmp == UNFEASIBLE_MIP_PROBLEM:
            return { 'status':'unfeasible' }
        elif tmp == UNBOUNDED_MIP_PROBLEM:
            return { 'status':'unbounded' }
        else:
            return { 'status':'optimized' }

    def optimizing_point(self):
        """
        Returns an optimal point for the MIP_Problem, if it exists.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> m = MIP_Problem()
        >>> m.add_space_dimensions_and_embed(2)
        >>> m.add_constraint(x >= 0)
        >>> m.add_constraint(y >= 0)
        >>> m.add_constraint(3 * x + 5 * y <= 10)
        >>> m.set_objective_function(x + y)
        >>> m.optimizing_point()
        point(10/3, 0/3)
        """
        cdef PPL_Generator *g
        sig_on()
        try:
            g = new_MIP_optimizing_point(self.thisptr[0])
        finally:
            sig_off()
        return _wrap_Generator(g[0])

    def OK(self):
        """
        Check if all the invariants are satisfied.

        OUTPUT:

        ``True`` if and only if ``self`` satisfies all the invariants.

        Examples:

        >>> from ppl import Variable, Constraint_System, MIP_Problem
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> m = MIP_Problem()
        >>> m.add_space_dimensions_and_embed(2)
        >>> m.add_constraint(x >= 0)
        >>> m.OK()
        True
        """
        return self.thisptr.OK()

cdef class MIP_Problem_constraints_iterator(object):
    """
    Wrapper for PPL's ``Constraint_System::const_iterator`` class.
    """
    def __cinit__(self, MIP_Problem pb):
        """
        The Cython constructor.

        See :class:`Constraint_System_iterator` for documentation.

        Tests:

        >>> from ppl import Constraint_System, Constraint_System_iterator
        >>> iter = Constraint_System_iterator( Constraint_System() )   # indirect doctest
        """
        self.pb = pb
        self.mip_csi_ptr = init_mip_cs_iterator(pb.thisptr[0])

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        delete_mip_cs_iterator(self.mip_csi_ptr)

    def __next__(Constraint_System_iterator self):
        r"""
        The next iteration.

        OUTPUT:

        A :class:`Generator`.

        Examples:

        >>> from ppl import Constraint_System, Variable, Constraint_System_iterator
        >>> x = Variable(0)
        >>> cs = Constraint_System( 5*x > 0 )
        >>> next(Constraint_System_iterator(cs))
        x0>0
        """
        if is_end_mip_cs_iterator((<MIP_Problem>self.pb).thisptr[0], self.mip_csi_ptr):
            raise StopIteration
        return _wrap_Constraint(next_mip_cs_iterator(self.mip_csi_ptr))

####################################################
### Polyhedron #####################################
####################################################
cdef class Polyhedron(_mutable_or_immutable):
    r"""
    Wrapper for PPL's ``Polyhedron`` class.

    An object of the class Polyhedron represents a convex polyhedron
    in the vector space.

    A polyhedron can be specified as either a finite system of
    constraints or a finite system of generators (see Section
    Representations of Convex Polyhedra) and it is always possible to
    obtain either representation. That is, if we know the system of
    constraints, we can obtain from this the system of generators that
    define the same polyhedron and vice versa. These systems can
    contain redundant members: in this case we say that they are not
    in the minimal form.

    INPUT/OUTPUT:

    This is an abstract base for :class:`C_Polyhedron` and
    :class:`NNC_Polyhedron`. You cannot instantiate this class.
    """
    def __init__(self):
        r"""
        The Python constructor.

        See also :class:`C_Polyhedron` and
        :class:`NNC_Polyhedron`. You must not instantiate
        :class:`Polyhedron` objects.

        Tests:

        >>> from ppl import Polyhedron
        >>> Polyhedron()
        Traceback (most recent call last):
        ...
        NotImplementedError: The Polyhedron class is abstract, you must not instantiate it.
        """
        raise NotImplementedError, 'The Polyhedron class is abstract, you must not instantiate it.'

    def __repr__(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        Examples:

        >>> from ppl import Variable, C_Polyhedron
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> C_Polyhedron( 5*x-2*y >=  x+y-1 )._repr_()
        'A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 ray, 1 line'

        Special cases::

        >>> C_Polyhedron(3, 'empty')._repr_()
        'The empty polyhedron in QQ^3'
        >>> C_Polyhedron(3, 'universe')._repr_()
        'The space-filling polyhedron in QQ^3'
        """
        dim = self.affine_dimension()
        ambient_dim = self.space_dimension()
        gs = self.minimized_generators()
        n_points = 0
        n_closure_points = 0
        n_lines = 0
        n_rays = 0
        for g in gs:
            if g.is_line():
                n_lines += 1
            elif g.is_ray():
                n_rays += 1
            elif g.is_point():
                n_points += 1
            elif g.is_closure_point():
                n_closure_points += 1
            else:
                raise RuntimeError
        if self.is_empty():
            return 'The empty polyhedron in QQ^'+str(ambient_dim)
        if self.is_universe():
            return 'The space-filling polyhedron in QQ^'+str(ambient_dim)
        desc = 'A ' + str(dim) + '-dimensional polyhedron'
        desc += ' in QQ'
        desc += '^' + str(ambient_dim)
        desc += ' defined as the convex hull of '
        first = True
        if n_points>0:
            if not first:
                desc += ", "
            first = False
            desc += str(n_points)
            if n_points==1:  desc += ' point'
            else:          desc += ' points'
        if n_closure_points>0:
            if not first:
                desc += ", "
            first = False
            desc += str(n_closure_points)
            if n_closure_points==1:  desc += ' closure_point'
            else:          desc += ' closure_points'
        if n_rays>0:
            if not first:
                desc += ", "
            first = False
            desc += str(n_rays)
            if n_rays==1:  desc += ' ray'
            else:          desc += ' rays'
        if n_lines>0:
            if not first:
                desc += ", "
            first = False
            desc += repr(n_lines)
            if n_lines==1: desc +=' line'
            else:          desc +=' lines'
        return desc;

    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUTPUT:

        Integer.

        Examples:

        >>> from ppl import Variable, C_Polyhedron
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> p = C_Polyhedron( 5*x-2*y >=  x+y-1 )
        >>> p.space_dimension()
        2
        """
        return self.thisptr.space_dimension()

    def affine_dimension(self):
        r"""
        Return the affine dimension of ``self``.

        OUTPUT:

        An integer. Returns 0 if ``self`` is empty. Otherwise, returns
        the affine dimension of ``self``.

        Examples:

        >>> from ppl import Variable, C_Polyhedron
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> p = C_Polyhedron( 5*x-2*y ==  x+y-1 )
        >>> p.affine_dimension()
        1
        """
        sig_on()
        cdef size_t dim = self.thisptr.affine_dimension()
        sig_off()
        return dim

    def constraints(self):
        r"""
        Returns the system of constraints.

        See also :meth:`minimized_constraints`.

        OUTPUT:

        A :class:`Constraint_System`.

        Examples:

        >>> from ppl import Variable, C_Polyhedron
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> p = C_Polyhedron(y >= 0)
        >>> p.add_constraint(x >= 0)
        >>> p.add_constraint(x+y >= 0)
        >>> p.constraints()
        Constraint_System {x1>=0, x0>=0, x0+x1>=0}
        >>> p.minimized_constraints()
        Constraint_System {x1>=0, x0>=0}
        """
        sig_on()
        cdef PPL_Constraint_System cs = self.thisptr.constraints()
        sig_off()
        return _wrap_Constraint_System(cs)

    def minimized_constraints(self):
        r"""
        Returns the minimized system of constraints.

        See also :meth:`constraints`.

        OUTPUT:

        A :class:`Constraint_System`.

        Examples:

        >>> from ppl import Variable, C_Polyhedron
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> p = C_Polyhedron(y >= 0)
        >>> p.add_constraint(x >= 0)
        >>> p.add_constraint(x+y >= 0)
        >>> p.constraints()
        Constraint_System {x1>=0, x0>=0, x0+x1>=0}
        >>> p.minimized_constraints()
        Constraint_System {x1>=0, x0>=0}
        """
        sig_on()
        cdef PPL_Constraint_System cs = self.thisptr.minimized_constraints()
        sig_off()
        return _wrap_Constraint_System(cs)

    def generators(self):
        r"""
        Returns the system of generators.

        See also :meth:`minimized_generators`.

        OUTPUT:

        A :class:`Generator_System`.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, point
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> p = C_Polyhedron(3,'empty')
        >>> p.add_generator(point(-x-y))
        >>> p.add_generator(point(0))
        >>> p.add_generator(point(+x+y))
        >>> p.generators()
        Generator_System {point(-1/1, -1/1, 0/1), point(0/1, 0/1, 0/1), point(1/1, 1/1, 0/1)}
        >>> p.minimized_generators()
        Generator_System {point(-1/1, -1/1, 0/1), point(1/1, 1/1, 0/1)}
        """
        sig_on()
        cdef PPL_Generator_System gs = self.thisptr.generators()
        sig_off()
        return _wrap_Generator_System(gs)

    def minimized_generators(self):
        r"""
        Returns the minimized system of generators.

        See also :meth:`generators`.

        OUTPUT:

        A :class:`Generator_System`.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, point
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> p = C_Polyhedron(3,'empty')
        >>> p.add_generator(point(-x-y))
        >>> p.add_generator(point(0))
        >>> p.add_generator(point(+x+y))
        >>> p.generators()
        Generator_System {point(-1/1, -1/1, 0/1), point(0/1, 0/1, 0/1), point(1/1, 1/1, 0/1)}
        >>> p.minimized_generators()
        Generator_System {point(-1/1, -1/1, 0/1), point(1/1, 1/1, 0/1)}
        """
        sig_on()
        cdef PPL_Generator_System gs = self.thisptr.minimized_generators()
        sig_off()
        return _wrap_Generator_System(gs)

    cdef _relation_with_generator(Polyhedron self, Generator g):
        r"""
        Helper method for :meth:`relation_with`.
        """
        rel = Poly_Gen_Relation(True)
        try:
            sig_on()
            try:
                rel.thisptr = new_relation_with(self.thisptr[0], g.thisptr[0])
            finally:
                sig_off()
        except BaseException:
            # rel.thisptr must be set to something valid or rel.__dealloc__() will segfault
            rel.thisptr = new PPL_Poly_Gen_Relation(PPL_Poly_Gen_Relation_nothing())
            raise
        return rel

    cdef _relation_with_constraint(Polyhedron self, Constraint c):
        r"""
        Helper method for :meth:`relation_with`.
        """
        rel = Poly_Con_Relation(True)
        try:
            sig_on()
            try:
                rel.thisptr = new_relation_with(self.thisptr[0], c.thisptr[0])
            finally:
                sig_off()
        except BaseException:
            # rel.thisptr must be set to something valid or rel.__dealloc__() will segfault
            rel.thisptr = new PPL_Poly_Con_Relation(PPL_Poly_Con_Relation_nothing())
            raise
        return rel

    def relation_with(self, arg):
        r"""
        Return the relations holding between the polyhedron ``self``
        and the generator or constraint ``arg``.

        INPUT:

        - ``arg`` -- a :class:`Generator` or a :class:`Constraint`.

        OUTPUT:

        A :class:`Poly_Gen_Relation` or a :class:`Poly_Con_Relation`
        according to the type of the input.

        Raises ``ValueError`` if ``self`` and the generator/constraint
        ``arg`` are dimension-incompatible.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, point, ray, Poly_Con_Relation
        >>> x = Variable(0);  y = Variable(1)
        >>> p = C_Polyhedron(2, 'empty')
        >>> p.add_generator( point(1*x+0*y) )
        >>> p.add_generator( point(0*x+1*y) )
        >>> p.minimized_constraints()
        Constraint_System {x0+x1-1==0, -x1+1>=0, x1>=0}
        >>> p.relation_with( point(1*x+1*y) )
        nothing
        >>> p.relation_with( point(1*x+1*y, 2) )
        subsumes
        >>> p.relation_with( x+y==-1 )
        is_disjoint
        >>> p.relation_with( x==y )
        strictly_intersects
        >>> p.relation_with( x+y<=1 )
        is_included, saturates
        >>> p.relation_with( x+y<1 )
        is_disjoint, saturates

        In a Python program you will usually use :meth:`relation_with`
        together with :meth:`~ppl.Poly_Gen_Relation.implies`
        or :meth:`~ppl.Poly_Con_Relation.implies`, for
        example:

        >>> p.relation_with( x+y<1 ).implies(Poly_Con_Relation.saturates())
        True

        You can only get relations with dimension-compatible
        generators or constraints:

        >>> z = Variable(2)
        >>> p.relation_with( point(x+y+z) )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::relation_with(g):
        this->space_dimension() == 2, g.space_dimension() == 3.
        >>> p.relation_with( z>0 )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::relation_with(c):
        this->space_dimension() == 2, c.space_dimension() == 3.
        """
        if isinstance(arg, Generator):
            return self._relation_with_generator(arg)
        if isinstance(arg, Constraint):
            return self._relation_with_constraint(arg)
        else:
            raise TypeError('Argument must be Generator or a Constraint')

    def is_empty(self):
        """
        Test if ``self`` is an empty polyhedron.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import C_Polyhedron
        >>> C_Polyhedron(3, 'empty').is_empty()
        True
        >>> C_Polyhedron(3, 'universe').is_empty()
        False
        """
        sig_on()
        cdef bint result = self.thisptr.is_empty()
        sig_off()
        return result

    def is_universe(self):
        """
        Test if ``self`` is a universe (space-filling) polyhedron.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import C_Polyhedron
        >>> C_Polyhedron(3, 'empty').is_universe()
        False
        >>> C_Polyhedron(3, 'universe').is_universe()
        True
        """
        sig_on()
        cdef bint result = self.thisptr.is_universe()
        sig_off()
        return result

    def is_topologically_closed(self):
        """
        Tests if ``self`` is topologically closed.

        OUTPUT:

        Returns ``True`` if and only if ``self`` is a topologically
        closed subset of the ambient vector space.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, NNC_Polyhedron
        >>> x = Variable(0);  y = Variable(1)
        >>> C_Polyhedron(3, 'universe').is_topologically_closed()
        True
        >>> C_Polyhedron( x>=1 ).is_topologically_closed()
        True
        >>> NNC_Polyhedron( x>1 ).is_topologically_closed()
        False
        """
        sig_on()
        cdef bint result = self.thisptr.is_topologically_closed()
        sig_off()
        return result

    def is_disjoint_from(self, Polyhedron y):
        r"""
        Tests whether ``self`` and ``y`` are disjoint.

        INPUT:

        - ``y`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` and ``y``
        are disjoint.

        Rayises a ``ValueError`` if ``self`` and ``y`` are
        topology-incompatible or dimension-incompatible.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, NNC_Polyhedron
        >>> x = Variable(0);  y = Variable(1)
        >>> C_Polyhedron(x<=0).is_disjoint_from( C_Polyhedron(x>=1) )
        True

        This is not allowed:

        >>> x = Variable(0);  y = Variable(1)
        >>> poly_1d = C_Polyhedron(x<=0)
        >>> poly_2d = C_Polyhedron(x+0*y>=1)
        >>> poly_1d.is_disjoint_from(poly_2d)
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::intersection_assign(y):
        this->space_dimension() == 1, y.space_dimension() == 2.

        Nor is this:

        >>> x = Variable(0);  y = Variable(1)
        >>> c_poly   =   C_Polyhedron( x<=0 )
        >>> nnc_poly = NNC_Polyhedron( x >0 )
        >>> c_poly.is_disjoint_from(nnc_poly)
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::intersection_assign(y):
        y is a NNC_Polyhedron.
        >>> NNC_Polyhedron(c_poly).is_disjoint_from(nnc_poly)
        True
        """
        cdef bint result
        sig_on()
        try:
            result = self.thisptr.is_disjoint_from(y.thisptr[0])
        finally:
            sig_off()
        return result

    def is_discrete(self):
        r"""
        Test whether ``self`` is discrete.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is discrete.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, point, ray
        >>> x = Variable(0);  y = Variable(1)
        >>> p = C_Polyhedron( point(1*x+2*y) )
        >>> p.is_discrete()
        True
        >>> p.add_generator( point(x) )
        >>> p.is_discrete()
        False
        """
        sig_on()
        cdef bint result = self.thisptr.is_discrete()
        sig_off()
        return result

    def is_bounded(self):
        r"""
        Test whether ``self`` is bounded.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is a bounded polyhedron.

        Examples:

        >>> from ppl import Variable, NNC_Polyhedron, point, closure_point, ray
        >>> x = Variable(0)
        >>> p = NNC_Polyhedron( point(0*x) )
        >>> p.add_generator( closure_point(1*x) )
        >>> p.is_bounded()
        True
        >>> p.add_generator( ray(1*x) )
        >>> p.is_bounded()
        False
        """
        sig_on()
        cdef bint result = self.thisptr.is_bounded()
        sig_off()
        return result

    def contains_integer_point(self):
        r"""
        Test whether ``self`` contains an integer point.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` contains an
        integer point.

        Examples:

        >>> from ppl import Variable, NNC_Polyhedron
        >>> x = Variable(0)
        >>> p = NNC_Polyhedron(x>0)
        >>> p.add_constraint(x<1)
        >>> p.contains_integer_point()
        False
        >>> p.topological_closure_assign()
        >>> p.contains_integer_point()
        True
        """
        sig_on()
        cdef bint result = self.thisptr.contains_integer_point()
        sig_off()
        return result

    def constrains(self, Variable var):
        r"""
        Test whether ``var`` is constrained in ``self``.

        INPUT:

        - ``var`` -- a :class:`Variable`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``var`` is
        constrained in ``self``.

        Raises a ``ValueError`` if ``var`` is not a space dimension of
        ``self``.

        Examples:

        >>> from ppl import Variable, C_Polyhedron
        >>> x = Variable(0)
        >>> p = C_Polyhedron(1, 'universe')
        >>> p.constrains(x)
        False
        >>> p = C_Polyhedron(x>=0)
        >>> p.constrains(x)
        True
        >>> y = Variable(1)
        >>> p.constrains(y)
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::constrains(v):
        this->space_dimension() == 1, v.space_dimension() == 2.
        """
        cdef bint result
        sig_on()
        try:
            result = self.thisptr.constrains(var.thisptr[0])
        finally:
            sig_off()
        return result

    def bounds_from_above(self, Linear_Expression expr):
        r"""
        Test whether the ``expr`` is bounded from above.

        INPUT:

        - ``expr`` -- a :class:`Linear_Expression`

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``expr`` is bounded
        from above in ``self``.

        Raises a ``ValueError`` if ``expr`` and ``this`` are
        dimension-incompatible.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, Linear_Expression
        >>> x = Variable(0);  y = Variable(1)
        >>> p = C_Polyhedron(y<=0)
        >>> p.bounds_from_above(x+1)
        False
        >>> p.bounds_from_above(Linear_Expression(y))
        True
        >>> p = C_Polyhedron(x<=0)
        >>> p.bounds_from_above(y+1)
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::bounds_from_above(e):
        this->space_dimension() == 1, e.space_dimension() == 2.
        """
        cdef bint result
        sig_on()
        try:
            result = self.thisptr.bounds_from_above(expr.thisptr[0])
        finally:
            sig_off()
        return result

    def bounds_from_below(self, Linear_Expression expr):
        r"""
        Test whether the ``expr`` is bounded from above.

        INPUT:

        - ``expr`` -- a :class:`Linear_Expression`

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``expr`` is bounded
        from above in ``self``.

        Raises a ``ValueError`` if ``expr`` and ``this`` are
        dimension-incompatible.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, Linear_Expression
        >>> x = Variable(0);  y = Variable(1)
        >>> p = C_Polyhedron(y>=0)
        >>> p.bounds_from_below(x+1)
        False
        >>> p.bounds_from_below(Linear_Expression(y))
        True
        >>> p = C_Polyhedron(x<=0)
        >>> p.bounds_from_below(y+1)
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::bounds_from_below(e):
        this->space_dimension() == 1, e.space_dimension() == 2.
        """
        cdef bint result
        sig_on()
        try:
            result = self.thisptr.bounds_from_below(expr.thisptr[0])
        finally:
            sig_off()
        return result

    def maximize(self, Linear_Expression expr):
        r"""
        Maximize ``expr``.

        INPUT:

        - ``expr`` -- a :class:`Linear_Expression`.

        OUTPUT:

        A dictionary with the following keyword:value pair:

        * ``'bounded'``: Boolean. Whether the linear expression
          ``expr`` is bounded from above on ``self``.

        If ``expr`` is bounded from above, the following additional
        keyword:value pairs are set to provide information about the
        supremum:

        * ``'sup_n'``: Integer. The numerator of the supremum value.

        * ``'sup_d'``: Non-zero integer. The denominator of the supremum
          value.

        * ``'maximum'``: Boolean. ``True`` if and only if the supremum
          is also the maximum value.

        * ``'generator'``: a :class:`Generator`. A point or closure
          point where expr reaches its supremum value.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, NNC_Polyhedron, Constraint_System, Linear_Expression
        >>> x = Variable(0);  y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert(x >= 0)
        >>> cs.insert(y >= 0)
        >>> cs.insert(3*x+5*y <= 10)
        >>> p = C_Polyhedron(cs)
        >>> pm = p.maximize(x+y)
        >>> for key in sorted(pm):
        ...     print key, pm[key]
        bounded True
        generator point(10/3, 0/3)
        maximum True
        sup_d 3
        sup_n 10

        Unbounded case::

        >>> cs = Constraint_System()
        >>> cs.insert(x > 0)
        >>> p = NNC_Polyhedron(cs)
        >>> p.maximize(+x)
        {'bounded': False}
        >>> pm = p.maximize(-x)
        >>> for key in pm:
        ...     print key, pm[key]
        bounded True
        generator closure_point(0/1)
        maximum False
        sup_d 1
        sup_n 0
        """
        cdef PPL_Coefficient sup_n
        cdef PPL_Coefficient sup_d
        cdef Generator g = point()
        cdef cppbool maximum
        sig_on()
        rc = self.thisptr.maximize(<PPL_Linear_Expression&>expr.thisptr[0], sup_n, sup_d, maximum, g.thisptr[0])
        sig_off()

        cdef long Int_sup_n = mpz_get_pyintlong(sup_n.get_mpz_t())
        cdef long Int_sup_d = mpz_get_pyintlong(sup_d.get_mpz_t())

        if rc:
            return { 'bounded':True, 'sup_n':Int_sup_n, 'sup_d':Int_sup_d, 'maximum':maximum, 'generator':g }
        else:
            return { 'bounded':False }

    def minimize(self, Linear_Expression expr):
        r"""
        Minimize ``expr``.

        INPUT:

        - ``expr`` -- a :class:`Linear_Expression`.

        OUTPUT:

        A dictionary with the following keyword:value pair:

        * ``'bounded'``: Boolean. Whether the linear expression
          ``expr`` is bounded from below on ``self``.

        If ``expr`` is bounded from below, the following additional
        keyword:value pairs are set to provide information about the
        infimum:

        * ``'inf_n'``: Integer. The numerator of the infimum value.

        * ``'inf_d'``: Non-zero integer. The denominator of the infimum
          value.

        * ``'minimum'``: Boolean. ``True`` if and only if the infimum
          is also the minimum value.

        * ``'generator'``: a :class:`Generator`. A point or closure
          point where expr reaches its infimum value.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, NNC_Polyhedron, Constraint_System, Linear_Expression
        >>> x = Variable(0);  y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x>=0 )
        >>> cs.insert( y>=0 )
        >>> cs.insert( 3*x+5*y<=10 )
        >>> p = C_Polyhedron(cs)
        >>> pm = p.minimize( x+y )
        >>> for key in sorted(pm):
        ...     print key, pm[key]
        bounded True
        generator point(0/1, 0/1)
        inf_d 1
        inf_n 0
        minimum True

        Unbounded case:

        >>> cs = Constraint_System()
        >>> cs.insert(x > 0)
        >>> p = NNC_Polyhedron(cs)
        >>> pm = p.minimize(+x)
        >>> for key in sorted(pm):
        ...    print key, pm[key]
        bounded True
        generator closure_point(0/1)
        inf_d 1
        inf_n 0
        minimum False
        >>> p.minimize( -x )
        {'bounded': False}
        """
        cdef PPL_Coefficient inf_n
        cdef PPL_Coefficient inf_d
        cdef Generator g = point()
        cdef cppbool minimum
        sig_on()
        rc = self.thisptr.minimize(<PPL_Linear_Expression&>expr.thisptr[0], inf_n, inf_d, minimum, g.thisptr[0])
        sig_off()

        cdef long Int_inf_n = mpz_get_pyintlong(inf_n.get_mpz_t())
        cdef long Int_inf_d = mpz_get_pyintlong(inf_d.get_mpz_t())

        if rc:
            return { 'bounded':True, 'inf_n':Int_inf_n, 'inf_d':Int_inf_d, 'minimum':minimum, 'generator':g }
        else:
            return { 'bounded':False }

    def contains(self, Polyhedron y):
        r"""
        Test whether ``self`` contains ``y``.

        INPUT:

        - ``y`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` contains ``y``.

        Raises a ``ValueError`` if ``self`` and ``y`` are
        topology-incompatible or dimension-incompatible.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, NNC_Polyhedron
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> p0 = C_Polyhedron( x>=0 )
        >>> p1 = C_Polyhedron( x>=1 )
        >>> p0.contains(p1)
        True
        >>> p1.contains(p0)
        False

        Errors are raised if the dimension or topology is not compatible:

        >>> p0.contains(C_Polyhedron(y>=0))
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::contains(y):
        this->space_dimension() == 1, y.space_dimension() == 2.
        >>> p0.contains(NNC_Polyhedron(x>0))
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::contains(y):
        y is a NNC_Polyhedron.
        """
        cdef bint result
        sig_on()
        try:
            result = self.thisptr.contains(y.thisptr[0])
        finally:
            sig_off()
        return result

    def strictly_contains(self, Polyhedron y):
        r"""
        Test whether ``self`` strictly contains ``y``.

        INPUT:

        - ``y`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` contains
        ``y`` and ``self`` does not equal ``y``.

        Raises a ``ValueError`` if ``self`` and ``y`` are
        topology-incompatible or dimension-incompatible.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, NNC_Polyhedron
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> p0 = C_Polyhedron( x>=0 )
        >>> p1 = C_Polyhedron( x>=1 )
        >>> p0.strictly_contains(p1)
        True
        >>> p1.strictly_contains(p0)
        False

        Errors are raised if the dimension or topology is not compatible:

        >>> p0.strictly_contains(C_Polyhedron(y>=0))
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::contains(y):
        this->space_dimension() == 1, y.space_dimension() == 2.
        >>> p0.strictly_contains(NNC_Polyhedron(x>0))
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::contains(y):
        y is a NNC_Polyhedron.
        """
        cdef bint result
        sig_on()
        try:
            result = self.thisptr.strictly_contains(y.thisptr[0])
        finally:
            sig_off()
        return result

    def add_constraint(self, Constraint c):
        r"""
        Add a constraint to the polyhedron.

        Adds a copy of constraint ``c`` to the system of constraints
        of ``self``, without minimizing the result.

        See alse :meth:`add_constraints`.

        INPUT:

        - ``c`` -- the :class:`Constraint` that will be added to the
          system of constraints of ``self``.

        OUTPUT:

        This method modifies the polyhedron ``self`` and does not
        return anything.

        Raises a ``ValueError`` if ``self`` and the constraint ``c`` are
        topology-incompatible or dimension-incompatible.

        Examples:

        >>> from ppl import Variable, C_Polyhedron
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> p = C_Polyhedron( y>=0 )
        >>> p.add_constraint( x>=0 )

         We just added a 1-d constraint to a 2-d polyhedron, this is
         fine. The other way is not:

        >>> p = C_Polyhedron( x>=0 )
        >>> p.add_constraint( y>=0 )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::add_constraint(c):
        this->space_dimension() == 1, c.space_dimension() == 2.

         The constraint must also be topology-compatible, that is,
         :class:`C_Polyhedron` only allows non-strict inequalities:

        >>> p = C_Polyhedron( x>=0 )
        >>> p.add_constraint( x< 1 )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::add_constraint(c):
        c is a strict inequality.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.add_constraint(c.thisptr[0])
        finally:
            sig_off()

    def add_generator(self, Generator g):
        r"""
        Add a generator to the polyhedron.

        Adds a copy of constraint ``c`` to the system of generators
        of ``self``, without minimizing the result.

        INPUT:

        - ``g`` -- the :class:`Generator` that will be added to the
          system of Generators of ``self``.

        OUTPUT:

        This method modifies the polyhedron ``self`` and does not
        return anything.

        Raises a ``ValueError`` if ``self`` and the generator ``g``
        are topology-incompatible or dimension-incompatible, or if
        ``self`` is an empty polyhedron and ``g`` is not a point.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, point, closure_point, ray
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> p = C_Polyhedron(1, 'empty')
        >>> p.add_generator( point(0*x) )

         We just added a 1-d generator to a 2-d polyhedron, this is
         fine. The other way is not:

        >>> p = C_Polyhedron(1, 'empty')
        >>> p.add_generator(  point(0*y) )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::add_generator(g):
        this->space_dimension() == 1, g.space_dimension() == 2.

         The constraint must also be topology-compatible, that is,
         :class:`C_Polyhedron` does not allow :func:`closure_point`
         generators:

        >>> p = C_Polyhedron( point(0*x+0*y) )
        >>> p.add_generator( closure_point(0*x) )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::add_generator(g):
        g is a closure point.

        Finally, ever non-empty polyhedron must have at least one
        point generator:

        >>> p = C_Polyhedron(3, 'empty')
        >>> p.add_generator( ray(x) )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::add_generator(g):
        *this is an empty polyhedron and g is not a point.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.add_generator(g.thisptr[0])
        finally:
            sig_off()

    def add_constraints(self, Constraint_System cs):
        r"""
        Add constraints to the polyhedron.

        Adds a copy of constraints in ``cs`` to the system of constraints
        of ``self``, without minimizing the result.

        See alse :meth:`add_constraint`.

        INPUT:

        - ``cs`` -- the :class:`Constraint_System` that will be added
          to the system of constraints of ``self``.

        OUTPUT:

        This method modifies the polyhedron ``self`` and does not
        return anything.

        Raises a ``ValueError`` if ``self`` and the constraints in
        ``cs`` are topology-incompatible or dimension-incompatible.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, Constraint_System
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert(x>=0)
        >>> cs.insert(y>=0)
        >>> p = C_Polyhedron( y<=1 )
        >>> p.add_constraints(cs)

        We just added a 1-d constraint to a 2-d polyhedron, this is
        fine. The other way is not:

        >>> p = C_Polyhedron( x<=1 )
        >>> p.add_constraints(cs)
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::add_recycled_constraints(cs):
        this->space_dimension() == 1, cs.space_dimension() == 2.

        The constraints must also be topology-compatible, that is,
        :class:`C_Polyhedron` only allows non-strict inequalities:

        >>> p = C_Polyhedron( x>=0 )
        >>> p.add_constraints( Constraint_System(x<0) )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::add_recycled_constraints(cs):
        cs contains strict inequalities.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.add_constraints(cs.thisptr[0])
        finally:
            sig_off()

    def add_generators(self, Generator_System gs):
        r"""
        Add generators to the polyhedron.

        Adds a copy of the generators in ``gs`` to the system of
        generators of ``self``, without minimizing the result.

        See alse :meth:`add_generator`.

        INPUT:

        - ``gs`` -- the :class:`Generator_System` that will be added
          to the system of constraints of ``self``.

        OUTPUT:

        This method modifies the polyhedron ``self`` and does not
        return anything.

        Raises a ``ValueError`` if ``self`` and one of the generators
        in ``gs`` are topology-incompatible or dimension-incompatible,
        or if ``self`` is an empty polyhedron and ``gs`` does not
        contain a point.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, Generator_System, point, ray, closure_point
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> gs = Generator_System()
        >>> gs.insert(point(0*x+0*y))
        >>> gs.insert(point(1*x+1*y))
        >>> p = C_Polyhedron(2, 'empty')
        >>> p.add_generators(gs)

        We just added a 1-d constraint to a 2-d polyhedron, this is
        fine. The other way is not:

        >>> p = C_Polyhedron(1, 'empty')
        >>> p.add_generators(gs)
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::add_recycled_generators(gs):
        this->space_dimension() == 1, gs.space_dimension() == 2.

        The constraints must also be topology-compatible, that is,
        :class:`C_Polyhedron` does not allow :func:`closure_point`
        generators:

        >>> p = C_Polyhedron( point(0*x+0*y) )
        >>> p.add_generators( Generator_System(closure_point(x) ))
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::add_recycled_generators(gs):
        gs contains closure points.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.add_generators(gs.thisptr[0])
        finally:
            sig_off()

    def unconstrain(self, Variable var):
        r"""
        Compute the cylindrification of ``self`` with respect to space
        dimension ``var``.

        INPUT:

        - ``var`` -- a :class:`Variable`. The space dimension that
          will be unconstrained.  Exceptions:

        OUTPUT:

        This method assigns the cylindrification to ``self`` and does
        not return anything.

        Raises a ``ValueError`` if ``var`` is not a space dimension of
        ``self``.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, point
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> p = C_Polyhedron( point(x+y) ); p
        A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point
        >>> p.unconstrain(x); p
        A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 line
        >>> z = Variable(2)
        >>> p.unconstrain(z)
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::unconstrain(var):
        this->space_dimension() == 2, required space dimension == 3.
        """
        sig_on()
        try:
            self.thisptr.unconstrain(var.thisptr[0])
        finally:
            sig_off()

    def intersection_assign(self, Polyhedron y):
        r"""
        Assign to ``self`` the intersection of ``self`` and ``y``.

        INPUT:

        - ``y`` -- a :class:`Polyhedron`

        OUTPUT:

        This method assigns the intersection to ``self`` and does not
        return anything.

        Raises a ``ValueError`` if ``self`` and and ``y`` are
        topology-incompatible or dimension-incompatible.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, NNC_Polyhedron
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> p = C_Polyhedron( 1*x+0*y >= 0 )
        >>> p.intersection_assign( C_Polyhedron(y>=0) )
        >>> p.constraints()
        Constraint_System {x0>=0, x1>=0}
        >>> z = Variable(2)
        >>> p.intersection_assign( C_Polyhedron(z>=0) )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::intersection_assign(y):
        this->space_dimension() == 2, y.space_dimension() == 3.
        >>> p.intersection_assign( NNC_Polyhedron(x+y<1) )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::intersection_assign(y):
        y is a NNC_Polyhedron.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.intersection_assign(y.thisptr[0])
        finally:
            sig_off()

    def poly_hull_assign(self, Polyhedron y):
        r"""
        Assign to ``self`` the poly-hull of ``self`` and ``y``.

        For any pair of NNC polyhedra `P_1` and `P_2`, the convex
        polyhedral hull (or poly-hull) of is the smallest NNC
        polyhedron that includes both `P_1` and `P_2`. The poly-hull
        of any pair of closed polyhedra in is also closed.

        INPUT:

        - ``y`` -- a :class:`Polyhedron`

        OUTPUT:

        This method assigns the poly-hull to ``self`` and does not
        return anything.

        Raises a ``ValueError`` if ``self`` and and ``y`` are
        topology-incompatible or dimension-incompatible.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, point, NNC_Polyhedron
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> p = C_Polyhedron( point(1*x+0*y) )
        >>> p.poly_hull_assign(C_Polyhedron( point(0*x+1*y) ))
        >>> p.generators()
        Generator_System {point(0/1, 1/1), point(1/1, 0/1)}

        ``self`` and ``y`` must be dimension- and topology-compatible,
        or an exception is raised:

        >>> z = Variable(2)
        >>> p.poly_hull_assign( C_Polyhedron(z>=0) )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::poly_hull_assign(y):
        this->space_dimension() == 2, y.space_dimension() == 3.
        >>> p.poly_hull_assign( NNC_Polyhedron(x+y<1) )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::poly_hull_assign(y):
        y is a NNC_Polyhedron.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.poly_hull_assign(y.thisptr[0])
        finally:
            sig_off()

    upper_bound_assign = poly_hull_assign

    def poly_difference_assign(self, Polyhedron y):
        r"""
        Assign to ``self`` the poly-difference of ``self`` and ``y``.

        For any pair of NNC polyhedra `P_1` and `P_2` the convex
        polyhedral difference (or poly-difference) of `P_1` and `P_2`
        is defined as the smallest convex polyhedron containing the
        set-theoretic difference `P_1\setminus P_2` of `P_1` and
        `P_2`.

        In general, even if `P_1` and `P_2` are topologically closed
        polyhedra, their poly-difference may be a convex polyhedron
        that is not topologically closed. For this reason, when
        computing the poly-difference of two :class:`C_Polyhedron`,
        the library will enforce the topological closure of the
        result.

        INPUT:

        - ``y`` -- a :class:`Polyhedron`

        OUTPUT:

        This method assigns the poly-difference to ``self`` and does
        not return anything.

        Raises a ``ValueError`` if ``self`` and and ``y`` are
        topology-incompatible or dimension-incompatible.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, point, closure_point, NNC_Polyhedron
        >>> x = Variable(0)
        >>> p = NNC_Polyhedron( point(0*x) )
        >>> p.add_generator( point(1*x) )
        >>> p.poly_difference_assign(NNC_Polyhedron( point(0*x) ))
        >>> p.minimized_constraints()
        Constraint_System {-x0+1>=0, x0>0}

        The poly-difference of :class:`C_polyhedron` is really its closure:

        >>> p = C_Polyhedron( point(0*x) )
        >>> p.add_generator( point(1*x) )
        >>> p.poly_difference_assign(C_Polyhedron( point(0*x) ))
        >>> p.minimized_constraints()
        Constraint_System {x0>=0, -x0+1>=0}

        ``self`` and ``y`` must be dimension- and topology-compatible,
        or an exception is raised:

        >>> y = Variable(1)
        >>> p.poly_difference_assign( C_Polyhedron(y>=0) )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::poly_difference_assign(y):
        this->space_dimension() == 1, y.space_dimension() == 2.
        >>> p.poly_difference_assign( NNC_Polyhedron(x+y<1) )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::poly_difference_assign(y):
        y is a NNC_Polyhedron.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.poly_difference_assign(y.thisptr[0])
        finally:
            sig_off()

    difference_assign = poly_difference_assign

    def drop_some_non_integer_points(self):
        r"""
        Possibly tighten ``self`` by dropping some points with
        non-integer coordinates.

        The modified polyhedron satisfies:

        * it is (not necessarily strictly) contained in the original
          polyhedron.

        * integral vertices (generating points with integer
          coordinates) of the original polyhedron are not removed.

        .. note::

            The modified polyhedron is not neccessarily a lattice
            polyhedron; Some vertices will, in general, still be
            rational. Lattice points interior to the polyhedron may be
            lost in the process.

        Examples:

        >>> from ppl import Variable, NNC_Polyhedron, Constraint_System
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x>=0 )
        >>> cs.insert( y>=0 )
        >>> cs.insert( 3*x+2*y<5 )
        >>> p = NNC_Polyhedron(cs)
        >>> p.minimized_generators()
        Generator_System {point(0/1, 0/1), closure_point(0/2, 5/2), closure_point(5/3, 0/3)}
        >>> p.drop_some_non_integer_points()
        >>> p.minimized_generators()
        Generator_System {point(0/1, 0/1), point(0/1, 2/1), point(4/3, 0/3)}
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        self.thisptr.drop_some_non_integer_points()
        sig_off()

    def topological_closure_assign(self):
        r"""
        Assign to ``self`` its topological closure.

        Examples:

        >>> from ppl import Variable, NNC_Polyhedron
        >>> x = Variable(0)
        >>> p = NNC_Polyhedron(x>0)
        >>> p.is_topologically_closed()
        False
        >>> p.topological_closure_assign()
        >>> p.is_topologically_closed()
        True
        >>> p.minimized_constraints()
        Constraint_System {x0>=0}
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        self.thisptr.topological_closure_assign()
        sig_off()

    def add_space_dimensions_and_embed(self, m):
        r"""
        Add ``m`` new space dimensions and embed ``self`` in the new
        vector space.

        The new space dimensions will be those having the highest
        indexes in the new polyhedron, which is characterized by a
        system of constraints in which the variables running through
        the new dimensions are not constrained. For instance, when
        starting from the polyhedron `P` and adding a third space
        dimension, the result will be the polyhedron

        .. MATH::

            \Big\{
            (x,y,z)^T \in \RR^3
            \Big|
            (x,y)^T \in P
            \Big\}

        INPUT:

        - ``m`` -- integer.

        OUTPUT:

        This method assigns the embedded polyhedron to ``self`` and
        does not return anything.

        Raises a ``ValueError`` if adding ``m`` new space dimensions
        would cause the vector space to exceed dimension
        ``self.max_space_dimension()``.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, point
        >>> x = Variable(0)
        >>> p = C_Polyhedron( point(3*x) )
        >>> p.add_space_dimensions_and_embed(1)
        >>> p.minimized_generators()
        Generator_System {line(0, 1), point(3/1, 0/1)}
        >>> p.add_space_dimensions_and_embed( p.max_space_dimension() )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::add_space_dimensions_and_embed(m):
        adding m new space dimensions exceeds the maximum allowed space dimension.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        m = int(m)
        sig_on()
        try:
            self.thisptr.add_space_dimensions_and_embed(m)
        finally:
            sig_off()

    def add_space_dimensions_and_project(self, m):
        r"""
        Add ``m`` new space dimensions and embed ``self`` in the new
        vector space.

        The new space dimensions will be those having the highest
        indexes in the new polyhedron, which is characterized by a
        system of constraints in which the variables running through
        the new dimensions are all constrained to be equal to `0`.
        For instance, when starting from the polyhedron `P` and adding
        a third space dimension, the result will be the polyhedron

        .. MATH::

            \Big\{
            (x,y,0)^T \in \RR^3
            \Big|
            (x,y)^T \in P
            \Big\}

        INPUT:

        - ``m`` -- integer.

        OUTPUT:

        This method assigns the projected polyhedron to ``self`` and
        does not return anything.

        Raises a ``ValueError`` if adding ``m`` new space dimensions
        would cause the vector space to exceed dimension
        ``self.max_space_dimension()``.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, point
        >>> x = Variable(0)
        >>> p = C_Polyhedron( point(3*x) )
        >>> p.add_space_dimensions_and_project(1)
        >>> p.minimized_generators()
        Generator_System {point(3/1, 0/1)}
        >>> p.add_space_dimensions_and_project( p.max_space_dimension() )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::add_space_dimensions_and_project(m):
        adding m new space dimensions exceeds the maximum allowed space dimension.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        m = int(m)
        sig_on()
        try:
            self.thisptr.add_space_dimensions_and_project(m)
        finally:
            sig_off()

    def concatenate_assign(self, Polyhedron y):
        r"""
        Assign to ``self`` the concatenation of ``self`` and ``y``.

        This functions returns the Cartiesian product of ``self`` and
        ``y``.

        Viewing a polyhedron as a set of tuples (its points), it is
        sometimes useful to consider the set of tuples obtained by
        concatenating an ordered pair of polyhedra. Formally, the
        concatenation of the polyhedra `P` and `Q` (taken in this
        order) is the polyhedron such that

        .. MATH::

            R =
            \Big\{
            (x_0,\dots,x_{n-1},y_0,\dots,y_{m-1})^T \in \RR^{n+m}
            \Big|
            (x_0,\dots,x_{n-1})^T \in P
            ,~
            (y_0,\dots,y_{m-1})^T \in Q
            \Big\}

        Another way of seeing it is as follows: first embed polyhedron
        `P` into a vector space of dimension `n+m` and then add a
        suitably renamed-apart version of the constraints defining
        `Q`.

        INPUT:

        - ``m`` -- integer.

        OUTPUT:

        This method assigns the concatenated polyhedron to ``self`` and
        does not return anything.

        Raises a ``ValueError`` if ``self`` and ``y`` are
        topology-incompatible or if adding ``y.space_dimension()`` new
        space dimensions would cause the vector space to exceed
        dimension ``self.max_space_dimension()``.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, NNC_Polyhedron, point
        >>> x = Variable(0)
        >>> p1 = C_Polyhedron( point(1*x) )
        >>> p2 = C_Polyhedron( point(2*x) )
        >>> p1.concatenate_assign(p2)
        >>> p1.minimized_generators()
        Generator_System {point(1/1, 2/1)}

        The polyhedra must be topology-compatible and not exceed the
        maximum space dimension:

        >>> p1.concatenate_assign( NNC_Polyhedron(1, 'universe') )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::concatenate_assign(y):
        y is a NNC_Polyhedron.
        >>> p1.concatenate_assign( C_Polyhedron(p1.max_space_dimension(), 'empty') )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::concatenate_assign(y):
        concatenation exceeds the maximum allowed space dimension.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        sig_on()
        try:
            self.thisptr.concatenate_assign(y.thisptr[0])
        finally:
            sig_off()

    def remove_higher_space_dimensions(self, new_dimension):
        r"""
        Remove the higher dimensions of the vector space so that the
        resulting space will have dimension ``new_dimension``.

        OUTPUT:

        This method modifies ``self`` and does not return anything.

        Raises a ``ValueError`` if ``new_dimensions`` is greater than
        the space dimension of ``self``.

        Examples:

        >>> from ppl import C_Polyhedron, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> p = C_Polyhedron(3*x+0*y==2)
        >>> p.remove_higher_space_dimensions(1)
        >>> p.minimized_constraints()
        Constraint_System {3*x0-2==0}
        >>> p.remove_higher_space_dimensions(2)
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::remove_higher_space_dimensions(nd):
        this->space_dimension() == 1, required space dimension == 2.
        """
        self.assert_mutable('The Polyhedron is not mutable!')
        new_dimension = int(new_dimension)
        sig_on()
        try:
            self.thisptr.remove_higher_space_dimensions(new_dimension)
        finally:
            sig_off()

    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        Examples:

        >>> cmd  = 'from ppl import C_Polyhedron, Variable\n'
        >>> cmd += 'x = Variable(0)\n'
        >>> cmd += 'y = Variable(1)\n'
        >>> cmd += 'p = C_Polyhedron(3*x+2*y==1)\n'
        >>> cmd += 'p.minimized_generators()\n'
        >>> cmd += 'p.ascii_dump()\n'
        >>> import subprocess
        >>> proc = subprocess.Popen(['python', '-c', cmd], stderr=subprocess.PIPE)
        >>> out, err = proc.communicate()
        >>> err
        'space_dim 2\n-ZE -EM  +CM +GM  +CS +GS  -CP -GP  -SC +SG \ncon_sys (up-to-date)\ntopology NECESSARILY_CLOSED\n2 x 2 SPARSE (sorted)\nindex_first_pending 2\nsize 3 -1 3 2 = (C)\nsize 3 1 0 0 >= (C)\n\ngen_sys (up-to-date)\ntopology NECESSARILY_CLOSED\n2 x 2 DENSE (not_sorted)\nindex_first_pending 2\nsize 3 0 2 -3 L (C)\nsize 3 2 0 1 P (C)\n\nsat_c\n0 x 0\n\nsat_g\n2 x 2\n0 0 \n0 1 \n\n'
        """
        sig_on()
        self.thisptr.ascii_dump()
        sig_off()

    def max_space_dimension(self):
        r"""
        Return the maximum space dimension all kinds of Polyhedron can handle.

        OUTPUT:

        Integer.

        Examples:

        >>> from ppl import C_Polyhedron
        >>> C_Polyhedron(1, 'empty').max_space_dimension()
        1152921504606846974
        """
        return self.thisptr.max_space_dimension()

    def OK(self, check_non_empty=False):
        """
        Check if all the invariants are satisfied.

        The check is performed so as to intrude as little as
        possible. If the library has been compiled with run-time
        assertions enabled, error messages are written on std::cerr in
        case invariants are violated. This is useful for the purpose
        of debugging the library.

        INPUT:

        - ``check_not_empty`` -- boolean. ``True`` if and only if, in
          addition to checking the invariants, ``self`` must be
          checked to be not empty.

        OUTPUT:

        ``True`` if and only if ``self`` satisfies all the invariants
        and either ``check_not_empty`` is ``False`` or ``self`` is not
        empty.

        Examples:

        >>> from ppl import Linear_Expression, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 3*x+2*y+1
        >>> e.OK()
        True
        """
        sig_on()
        cdef bint result = self.thisptr.OK()
        sig_off()
        return result

    def __hash__(self):
        r"""
        Hash value for polyhedra.

        Tests:

        >>> from ppl import Constraint_System, Variable, C_Polyhedron
        >>> x = Variable(0)
        >>> p = C_Polyhedron( 5*x >= 3 )
        >>> p.set_immutable()
        >>> hash(p)
        1

        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0 )
        >>> cs.insert( y >= 0 )
        >>> p = C_Polyhedron(cs)
        >>> p.set_immutable()
        >>> hash(p)
        2

        >>> hash(C_Polyhedron(x >= 0))
        Traceback (most recent call last):
        ...
        TypeError: mutable polyhedra are unhashable
        """
        if self.is_mutable():
            raise TypeError("mutable polyhedra are unhashable")
        # TODO: the hash code from PPL looks like being the dimension!
        return self.thisptr[0].hash_code()

    def __richcmp__(Polyhedron lhs, Polyhedron rhs, op):
        r"""
        Comparison for polyhedra.

        INPUT:

        - ``lhs``, ``rhs`` -- :class:`Polyhedron`.

        - ``op`` -- integer. The comparison operation to be performed.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import Variable, C_Polyhedron
        >>> x = Variable(0)
        >>> C_Polyhedron(x >= 0) > C_Polyhedron(x >= 1)
        True
        """
        cdef result
        sig_on()
        if op==0:      # <   0
            result = rhs.strictly_contains(lhs)
        elif op==1:    # <=  1
            result = rhs.contains(lhs)
        elif op==2:    # ==  2
            result = (lhs.thisptr[0] == rhs.thisptr[0])
        elif op==4:    # >   4
            result = lhs.strictly_contains(rhs)
        elif op==5:    # >=  5
            result = lhs.contains(rhs)
        elif op==3:    # !=  3
            result = (lhs.thisptr[0] != rhs.thisptr[0])
        else:
            raise RuntimeError  # unreachable
        sig_off()
        return result



####################################################
### C_Polyhedron ###################################
####################################################
cdef class C_Polyhedron(Polyhedron):
    r"""
    Wrapper for PPL's ``C_Polyhedron`` class.

    An object of the class :class:`C_Polyhedron` represents a
    topologically closed convex polyhedron in the vector space. See
    :class:`NNC_Polyhedron` for more general (not necessarily closed)
    polyhedra.

    When building a closed polyhedron starting from a system of
    constraints, an exception is thrown if the system contains a
    strict inequality constraint. Similarly, an exception is thrown
    when building a closed polyhedron starting from a system of
    generators containing a closure point.

    INPUT:

    - ``arg`` -- the defining data of the polyhedron. Any one of the
      following is accepted:

      * A non-negative integer. Depending on ``degenerate_element``,
        either the space-filling or the empty polytope in the given
        dimension ``arg`` is constructed.

      * A :class:`Constraint_System`.

      * A :class:`Generator_System`.

      * A single :class:`Constraint`.

      * A single :class:`Generator`.

      * A :class:`C_Polyhedron`.

    - ``degenerate_element`` -- string, either ``'universe'`` or
      ``'empty'``. Only used if ``arg`` is an integer.

    OUTPUT:

    A :class:`C_Polyhedron`.

    Examples:

    >>> from ppl import Constraint, Constraint_System, Generator, Generator_System, Variable, C_Polyhedron, point, ray
    >>> x = Variable(0)
    >>> y = Variable(1)
    >>> C_Polyhedron( 5*x-2*y >=  x+y-1 )
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 ray, 1 line
    >>> cs = Constraint_System()
    >>> cs.insert( x >= 0 )
    >>> cs.insert( y >= 0 )
    >>> C_Polyhedron(cs)
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 2 rays
    >>> C_Polyhedron( point(x+y) )
    A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point
    >>> gs = Generator_System()
    >>> gs.insert( point(-x-y) )
    >>> gs.insert( ray(x) )
    >>> C_Polyhedron(gs)
    A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 ray

    The empty and universe polyhedra are constructed like this:

    >>> C_Polyhedron(3, 'empty')
    The empty polyhedron in QQ^3
    >>> C_Polyhedron(3, 'empty').constraints()
    Constraint_System {-1==0}
    >>> C_Polyhedron(3, 'universe')
    The space-filling polyhedron in QQ^3
    >>> C_Polyhedron(3, 'universe').constraints()
    Constraint_System {}

    Note that, by convention, the generator system of a polyhedron is
    either empty or contains at least one point. In particular, if you
    define a polyhedron via a non-empty :class:`Generator_System` it
    must contain a point (at any position). If you start with a single
    generator, this generator must be a point:

    >>> C_Polyhedron( ray(x) )
    Traceback (most recent call last):
    ...
    ValueError: PPL::C_Polyhedron::C_Polyhedron(gs):
    *this is an empty polyhedron and
    the non-empty generator system gs contains no points.
    """
    def __cinit__(self, arg, degenerate_element='universe'):
        """
        The Cython constructor.

        See :class:`C_Polyhedron` for documentation.

        Tests:

        >>> from ppl import C_Polyhedron
        >>> C_Polyhedron(3, 'empty')   # indirect doctest
        The empty polyhedron in QQ^3
        """
        if isinstance(arg, C_Polyhedron):
            ph = <C_Polyhedron>arg
            self.thisptr = new PPL_C_Polyhedron(<PPL_C_Polyhedron&>ph.thisptr[0])
            return
        if isinstance(arg, Generator):
            arg = Generator_System(arg)
        if isinstance(arg, Constraint):
            arg = Constraint_System(arg)
        if isinstance(arg, Generator_System):
            gs = <Generator_System>arg
            self.thisptr = new PPL_C_Polyhedron(gs.thisptr[0])
            return
        if isinstance(arg, Constraint_System):
            cs = <Constraint_System>arg
            self.thisptr = new PPL_C_Polyhedron(cs.thisptr[0])
            return
        try:
            dim = int(arg)
            assert dim>=0
        except ValueError:
            raise ValueError, 'Cannot initialize C_Polyhedron with '+str(arg)+'.'
        degenerate_element = degenerate_element.lower()
        if degenerate_element=='universe':
            self.thisptr = new PPL_C_Polyhedron(<PPL_dimension_type>dim, UNIVERSE)
            return
        elif degenerate_element=='empty':
            self.thisptr = new PPL_C_Polyhedron(<PPL_dimension_type>dim, EMPTY)
            return
        else:
            raise ValueError, 'Unknown value: degenerate_element='+str(degenerate_element)+'.'

    def __init__(self, *args):
        """
        The Python destructor.

        See :class:`C_Polyhedron` for documentation.

        Tests:

        >>> from ppl import C_Polyhedron
        >>> C_Polyhedron(3, 'empty')   # indirect doctest
        The empty polyhedron in QQ^3
        """
        # override Polyhedron.__init__
        pass

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr

    def __reduce__(self):
        """
        Pickle object

        Tests:

        >>> from ppl import C_Polyhedron, Variable
        >>> P = C_Polyhedron(3, 'empty')
        >>> loads(dumps(P))
        The empty polyhedron in QQ^3

        >>> Q = C_Polyhedron(5, 'universe')
        >>> loads(dumps(Q))
        The space-filling polyhedron in QQ^5

        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> H = C_Polyhedron( 5*x-2*y >=  x+y-1 )
        >>> loads(dumps(H))
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 ray, 1 line
        """
        if self.is_empty():
            return (C_Polyhedron, (self.space_dimension(), 'empty'))
        elif self.is_universe():
            return (C_Polyhedron, (self.space_dimension(), 'universe'))
        else:
            return (C_Polyhedron, (self.generators(),))


cdef class NNC_Polyhedron(Polyhedron):
    r"""
    Wrapper for PPL's ``NNC_Polyhedron`` class.

    An object of the class ``NNC_Polyhedron`` represents a not
    necessarily closed (NNC) convex polyhedron in the vector space.

    Note: Since NNC polyhedra are a generalization of closed
    polyhedra, any object of the class :class:`C_Polyhedron` can be
    (explicitly) converted into an object of the class
    :class:`NNC_Polyhedron`. The reason for defining two different
    classes is that objects of the class :class:`C_Polyhedron` are
    characterized by a more efficient implementation, requiring less
    time and memory resources.

    INPUT:

    - ``arg`` -- the defining data of the polyhedron. Any one of the
      following is accepted:

      * An non-negative integer. Depending on ``degenerate_element``,
        either the space-filling or the empty polytope in the given
        dimension ``arg`` is constructed.

      * A :class:`Constraint_System`.

      * A :class:`Generator_System`.

      * A single :class:`Constraint`.

      * A single :class:`Generator`.

      * A :class:`NNC_Polyhedron`.

      * A :class:`C_Polyhedron`.

    - ``degenerate_element`` -- string, either ``'universe'`` or
      ``'empty'``. Only used if ``arg`` is an integer.

    OUTPUT:

    A :class:`C_Polyhedron`.

    Examples:

    >>> from ppl import Constraint, Constraint_System, Generator, Generator_System, Variable, NNC_Polyhedron, point, ray, closure_point
    >>> x = Variable(0)
    >>> y = Variable(1)
    >>> NNC_Polyhedron( 5*x-2*y >  x+y-1 )
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 closure_point, 1 ray, 1 line
    >>> cs = Constraint_System()
    >>> cs.insert( x > 0 )
    >>> cs.insert( y > 0 )
    >>> NNC_Polyhedron(cs)
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 closure_point, 2 rays
    >>> NNC_Polyhedron( point(x+y) )
    A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point
    >>> gs = Generator_System()
    >>> gs.insert( point(-y) )
    >>> gs.insert( closure_point(-x-y) )
    >>> gs.insert( ray(x) )
    >>> p = NNC_Polyhedron(gs); p
    A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 closure_point, 1 ray
    >>> p.minimized_constraints()
    Constraint_System {x1+1==0, x0+1>0}

    Note that, by convention, every polyhedron must contain a point:

    >>> NNC_Polyhedron( closure_point(x+y) )
    Traceback (most recent call last):
    ...
    ValueError: PPL::NNC_Polyhedron::NNC_Polyhedron(gs):
    *this is an empty polyhedron and
    the non-empty generator system gs contains no points.
    """
    def __cinit__(self, arg, degenerate_element='universe'):
        """
        The Cython constructor.

        See :class:`NNC_Polyhedron` for documentation.

        Tests:

        >>> from ppl import NNC_Polyhedron
        >>> NNC_Polyhedron(3, 'empty')   # indirect doctest
        The empty polyhedron in QQ^3
        """
        if isinstance(arg, NNC_Polyhedron):
            p_nnc = <NNC_Polyhedron>arg
            self.thisptr = new PPL_NNC_Polyhedron(<PPL_NNC_Polyhedron&>p_nnc.thisptr[0])
            return
        if isinstance(arg, C_Polyhedron):
            p_c = <C_Polyhedron>arg
            self.thisptr = new PPL_NNC_Polyhedron(<PPL_C_Polyhedron&>p_c.thisptr[0])
            return
        if isinstance(arg, Generator):
            arg = Generator_System(arg)
        if isinstance(arg, Constraint):
            arg = Constraint_System(arg)
        if isinstance(arg, Generator_System):
            gs = <Generator_System>arg
            self.thisptr = new PPL_NNC_Polyhedron(gs.thisptr[0])
            return
        if isinstance(arg, Constraint_System):
            cs = <Constraint_System>arg
            self.thisptr = new PPL_NNC_Polyhedron(cs.thisptr[0])
            return
        try:
            dim = int(arg)
            assert dim>=0
        except ValueError:
            raise ValueError('Cannot initialize NNC_Polyhedron with '+str(arg)+'.')
        degenerate_element = degenerate_element.lower()
        if degenerate_element=='universe':
            self.thisptr = new PPL_NNC_Polyhedron(<PPL_dimension_type>dim, UNIVERSE)
            return
        elif degenerate_element=='empty':
            self.thisptr = new PPL_NNC_Polyhedron(<PPL_dimension_type>dim, EMPTY)
            return
        else:
            raise ValueError('Unknown value: degenerate_element='+str(degenerate_element)+'.')

    def __init__(self, *args):
        """
        The Python destructor.

        See :class:`NNC_Polyhedron` for documentation.

        Tests:

        >>> from ppl import NNC_Polyhedron
        >>> NNC_Polyhedron(3, 'empty')   # indirect doctest
        The empty polyhedron in QQ^3
        """
        # override Polyhedron.__init__
        pass

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr

    def __reduce__(self):
        """
        Pickle object

        Tests:

        >>> from ppl import NNC_Polyhedron, Variable
        >>> P = NNC_Polyhedron(3, 'empty')
        >>> loads(dumps(P))
        The empty polyhedron in QQ^3

        >>> Q = NNC_Polyhedron(5, 'universe')
        >>> loads(dumps(Q))
        The space-filling polyhedron in QQ^5

        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> H = NNC_Polyhedron( 5*x-2*y >  x+y-1 )
        >>> loads(dumps(H))
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point,
        1 closure_point, 1 ray, 1 line
        """
        if self.is_empty():
            return (NNC_Polyhedron, (self.space_dimension(), 'empty'))
        elif self.is_universe():
            return (NNC_Polyhedron, (self.space_dimension(), 'universe'))
        else:
            return (NNC_Polyhedron, (self.generators(),))


cdef class Variable(object):
    r"""
    Wrapper for PPL's ``Variable`` class.

    A dimension of the vector space.

    An object of the class Variable represents a dimension of the space, that is
    one of the Cartesian axes. Variables are used as basic blocks in order to
    build more complex linear expressions. Each variable is identified by a
    non-negative integer, representing the index of the corresponding Cartesian
    axis (the first axis has index 0). The space dimension of a variable is the
    dimension of the vector space made by all the Cartesian axes having an index
    less than or equal to that of the considered variable; thus, if a variable
    has index `i`, its space dimension is `i+1`.

    INPUT:

    - ``i`` -- integer. The index of the axis.

    OUTPUT:

    A :class:`Variable`

    Examples:

    >>> from ppl import Variable
    >>> x = Variable(123)
    >>> x.id()
    123
    >>> x
    x123

    Note that the "meaning" of an object of the class Variable is completely
    specified by the integer index provided to its constructor: be careful not
    to be mislead by C++ language variable names. For instance, in the following
    example the linear expressions ``e1`` and ``e2`` are equivalent, since the
    two variables ``x`` and ``z`` denote the same Cartesian axis:

    >>> x = Variable(0)
    >>> y = Variable(1)
    >>> z = Variable(0)
    >>> e1 = x + y; e1
    x0+x1
    >>> e2 = y + z; e2
    x0+x1
    >>> e1 - e2
    0
    """
    def __cinit__(self, PPL_dimension_type i):
        """
        The Cython constructor.

        See :class:`Variable` for documentation.

        Tests:

        >>> from ppl import Variable
        >>> Variable(123)   # indirect doctest
        x123
        """
        self.thisptr = new PPL_Variable(i)

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr

    def id(self):
        """
        Return the index of the Cartesian axis associated to the variable.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(123)
        >>> x.id()
        123
        """
        return self.thisptr.id()

    def OK(self):
        """
        Checks if all the invariants are satisfied.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> x.OK()
        True
        """
        return self.thisptr.OK()

    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUPUT:

        Integer. The returned value is ``self.id()+1``.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> x.space_dimension()
        1
        """
        return self.thisptr.space_dimension()

    def __repr__(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> x.__repr__()
        'x0'
        """
        return 'x{0}'.format(self.id())

    def __add__(self, other):
        r"""
        Return the sum ``self`` + ``other``.

        INPUT:

        - ``self``, ``other`` -- anything convertible to
          ``Linear_Expression``: An integer, a :class:`Variable`, or a
          :class:`Linear_Expression`.

        OUTPUT:

        A :class:`Linear_Expression` representing ``self`` + ``other``.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0); y = Variable(1)
        >>> x + 15
        x0+15
        >>> 15 + y
        x1+15
        """
        return Linear_Expression(self)+Linear_Expression(other)

    def __sub__(self, other):
        r"""
        Return the difference ``self`` - ``other``.

        INPUT:

        - ``self``, ``other`` -- anything convertible to
          ``Linear_Expression``: An integer, a :class:`Variable`, or a
          :class:`Linear_Expression`.

        OUTPUT:

        A :class:`Linear_Expression` representing ``self`` - ``other``.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0); y = Variable(1)
        >>> x - 15
        x0-15
        >>> 15 - y
        -x1+15
        """
        return Linear_Expression(self)-Linear_Expression(other)

    def __mul__(self, other):
        r"""
        Return the product ``self`` * ``other``.

        INPUT:

        - ``self``, ``other`` -- One must be an integer, the other a
          :class:`Variable`.

        OUTPUT:

        A :class:`Linear_Expression` representing ``self`` * ``other``.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0); y = Variable(1)
        >>> x * 15
        15*x0
        >>> 15 * y
        15*x1
        """
        if isinstance(self, Variable):
            return Linear_Expression(self) * other
        else:
            return Linear_Expression(other) * self

    def __pos__(self):
        r"""
        Return ``self`` as :class:`Linear_Expression`

        OUTPUT:

        The :class:`Linear_Expression` ``+self``

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0); x
        x0
        >>> +x
        x0
        """
        return Linear_Expression(self)

    def __neg__(self):
        r"""
        Return -``self`` as :class:`Linear_Expression`

        OUTPUT:

        The :class:`Linear_Expression` ``-self``

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0); x
        x0
        >>> -x
        -x0
        """
        return Linear_Expression(self)*(-1)

    def __richcmp__(self, other, op):
        """
        Construct :class:`Constraint` from equalities or inequalities.

        INPUT:

        - ``self``, ``other`` -- anything convertible to a
          :class:`Linear_Expression`

        - ``op`` -- the operation.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> x <  y
        -x0+x1>0
        >>> x <= 0
        -x0>=0
        >>> x == y-y
        x0==0
        >>> x >= -2
        x0+2>=0
        >>> x >  0
        x0>0
        >>> 0 == 1    # watch out!
        False
        >>> 0*x == 1
        -1==0
        """
        return _make_Constraint_from_richcmp(self, other, op)


####################################################
### Linear_Expression ##############################
####################################################
cdef class Linear_Expression(object):
    r"""
    Wrapper for PPL's ``PPL_Linear_Expression`` class.

    INPUT:

    The constructor accepts zero, one, or two arguments.

    If there are two arguments ``Linear_Expression(a,b)``, they are
    interpreted as

    - ``a`` -- an iterable of integer coefficients, for example a
      list.

    - ``b`` -- an integer. The inhomogeneous term.

    A single argument ``Linear_Expression(arg)`` is interpreted as

    - ``arg`` -- something that determines a linear
      expression. Possibilities are:

      * a :class:`Variable`: The linear expression given by that
        variable.

      * a :class:`Linear_Expression`: The copy constructor.

      * an integer: Constructs the constant linear expression.

    No argument is the default constructor and returns the zero linear
    expression.

    OUTPUT:

    A :class:`Linear_Expression`

    Examples:

    >>> from ppl import Variable, Linear_Expression
    >>> Linear_Expression([1,2,3,4],5)
    x0+2*x1+3*x2+4*x3+5
    >>> Linear_Expression(10)
    10
    >>> Linear_Expression()
    0
    >>> Linear_Expression(10).inhomogeneous_term()
    10
    >>> x = Variable(123)
    >>> expr = x+1; expr
    x123+1
    >>> expr.OK()
    True
    >>> expr.coefficient(x)
    1
    >>> expr.coefficient( Variable(124) )
    0
    """
    def __cinit__(self, *args):
        """
        The Cython constructor.

        See :class:`Linear_Expression` for documentation.

        Tests:

        >>> from ppl import Linear_Expression
        >>> Linear_Expression(10)   # indirect doctest
        10
        """
        cdef long i
        if len(args)==2:
            a = args[0]
            b = args[1]
            ex = Linear_Expression(0)
            for i in range(len(a)):
                ex += Variable(i) * a[i]
            arg = ex + b
        elif len(args) == 1:
            arg = args[0]
        elif len(args) == 0:
            self.thisptr = new PPL_Linear_Expression()
            return
        else:
            raise ValueError("Cannot initialize with more than 2 arguments.")

        if isinstance(arg, Variable):
            v = <Variable>arg
            self.thisptr = new PPL_Linear_Expression(v.thisptr[0])
            return
        if isinstance(arg, Linear_Expression):
            e = <Linear_Expression>arg
            self.thisptr = new PPL_Linear_Expression(e.thisptr[0])
            return
        self.thisptr = new PPL_Linear_Expression(PPL_Coefficient_from_pyobject(arg))

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr

    def space_dimension(self):
        """
        Return the dimension of the vector space necessary for the
        linear expression.

        OUTPUT:

        Integer.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> (x+y+1).space_dimension()
        2
        >>> (x+y).space_dimension()
        2
        >>> (y+1).space_dimension()
        2
        >>> (x+1).space_dimension()
        1
        >>> (y+1-y).space_dimension()
        2
        """
        return self.thisptr.space_dimension()

    def coefficient(self, Variable v):
        """
        Return the coefficient of the variable ``v``.

        INPUT:

        - ``v`` -- a :class:`Variable`.

        OUTPUT:

        An integer.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> e = 3*x+1
        >>> e.coefficient(x)
        3
        """
        return mpz_get_pyintlong(self.thisptr.coefficient(v.thisptr[0]).get_mpz_t())

    def coefficients(self):
        """
        Return the coefficients of the linear expression.

        OUTPUT:

        A tuple of integers of length :meth:`space_dimension`.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0);  y = Variable(1)
        >>> e = 3*x+5*y+1
        >>> e.coefficients()
        (3, 5)
        """
        cdef int d = self.space_dimension()
        cdef int i
        cdef list coeffs = [None]*d
        for i in range(d):
            coeffs[i] = mpz_get_pyintlong(self.thisptr.coefficient(PPL_Variable(i)).get_mpz_t())
        return tuple(coeffs)

    def inhomogeneous_term(self):
        """
        Return the inhomogeneous term of the linear expression.

        OUTPUT:

        Integer.

        Examples:

        >>> from ppl import Variable, Linear_Expression
        >>> Linear_Expression(10).inhomogeneous_term()
        10
        """
        return mpz_get_pyintlong(self.thisptr.inhomogeneous_term().get_mpz_t())

    def is_zero(self):
        """
        Test if ``self`` is the zero linear expression.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import Variable, Linear_Expression
        >>> Linear_Expression(0).is_zero()
        True
        >>> Linear_Expression(10).is_zero()
        False
        """
        return self.thisptr.is_zero()

    def all_homogeneous_terms_are_zero(self):
        """
        Test if ``self`` is a constant linear expression.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import Variable, Linear_Expression
        >>> Linear_Expression(10).all_homogeneous_terms_are_zero()
        True
        """
        return self.thisptr.all_homogeneous_terms_are_zero()

    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        Examples:

        >>> cmd  = 'from ppl import Linear_Expression, Variable\n'
        >>> cmd += 'x = Variable(0)\n'
        >>> cmd += 'y = Variable(1)\n'
        >>> cmd += 'e = 3*x+2*y+1\n'
        >>> cmd += 'e.ascii_dump()\n'
        >>> import subprocess
        >>> proc = subprocess.Popen(['python', '-c', cmd], stderr=subprocess.PIPE)
        >>> out, err = proc.communicate()
        >>> print err
        size 3 1 3 2
        """
        self.thisptr.ascii_dump()

    def OK(self):
        """
        Check if all the invariants are satisfied.

        Examples:

        >>> from ppl import Linear_Expression, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 3*x+2*y+1
        >>> e.OK()
        True
        """
        return self.thisptr.OK()

    def __repr__(self):
        r"""
        Return a string representation of the linear expression.

        OUTPUT:

        A string.

        Examples:

        >>> from ppl import Linear_Expression, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> x+1
        x0+1
        >>> x+1-x
        1
        >>> 2*x
        2*x0
        >>> x-x-1
        -1
        >>> x-x
        0
        """
        s = ''
        first = True
        for i in range(0,self.space_dimension()):
            x = Variable(i)
            coeff = self.coefficient(x)
            if coeff==0: continue
            if first and coeff==1:
                s += '%r' % x
                first = False
            elif first and coeff==-1:
                s += '-%r' % x
                first = False
            elif first and coeff!=1:
                s += '%d*%r' % (coeff, x)
                first = False
            elif coeff==1:
                s += '+%r' % x
            elif coeff==-1:
                s += '-%r' % x
            else:
                s += '%+d*%r' % (coeff, x)
        inhomog = self.inhomogeneous_term()
        if inhomog != 0:
            if first:
                s += '%d' % inhomog
                first = False
            else:
                s += '%+d' % inhomog
        if first:
            s = '0'
        return s

    def __add__(self, other):
        r"""
        Add ``self`` and ``other``.

        INPUT:

        - ``self``, ``other`` -- anything that can be used to
          construct a :class:`Linear_Expression`. One of them, not
          necessarily ``self``, is guaranteed to be a
          :class:``Linear_Expression``, otherwise Python would not
          have called this method.

        OUTPUT:

        The sum as a :class:`Linear_Expression`

        Examples:

        >>> from ppl import Linear_Expression, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> 9+x+y+(1+x)+y+y
        2*x0+3*x1+10
        """
        cdef Linear_Expression lhs = Linear_Expression(self)
        cdef Linear_Expression rhs = Linear_Expression(other)
        cdef Linear_Expression result = Linear_Expression()
        result.thisptr[0] = lhs.thisptr[0] + rhs.thisptr[0]
        return result

    def __sub__(self, other):
        r"""
        Subtract ``other`` from ``self``.

        INPUT:

        - ``self``, ``other`` -- anything that can be used to
          construct a :class:`Linear_Expression`. One of them, not
          necessarily ``self``, is guaranteed to be a
          :class:``Linear_Expression``, otherwise Python would not
          have called this method.

        OUTPUT:

        The difference as a :class:`Linear_Expression`

        Examples:

        >>> from ppl import Linear_Expression, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> 9-x-y-(1-x)-y-y
        -3*x1+8
        """
        cdef Linear_Expression lhs = Linear_Expression(self)
        cdef Linear_Expression rhs = Linear_Expression(other)
        cdef Linear_Expression result = Linear_Expression()
        result.thisptr[0] = lhs.thisptr[0] - rhs.thisptr[0]
        return result

    def __mul__(self, other):
        r"""
        Multiply ``self`` with ``other``.

        INPUT:

        - ``self``, ``other`` -- anything that can be used to
          construct a :class:`Linear_Expression`. One of them, not
          necessarily ``self``, is guaranteed to be a
          :class:``Linear_Expression``, otherwise Python would not
          have called this method.

        OUTPUT:

        The product as a :class:`Linear_Expression`

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> 8*(x+1)
        8*x0+8
        >>> y*8
        8*x1
        >>> 2**128 * x
        340282366920938463463374607431768211456*x0
        """
        cdef Linear_Expression e
        cdef c
        if isinstance(self, Linear_Expression):
            e = <Linear_Expression>self
            c = other
        else:
            e = <Linear_Expression>other
            c = self

        cdef PPL_Coefficient cc = PPL_Coefficient_from_pyobject(c)
        cdef Linear_Expression result = Linear_Expression()
        result.thisptr[0] = e.thisptr[0] * cc
        return result

    def __pos__(self):
        """
        Return ``self``.

        Examples:

        >>> from ppl import Variable, Linear_Expression
        >>> +Linear_Expression(1)
        1
        >>> x = Variable(0)
        >>> +(x+1)
        x0+1
        """
        return self

    def __neg__(self):
        """
        Return the negative of ``self``.

        Examples:

        >>> from ppl import Variable, Linear_Expression
        >>> -Linear_Expression(1)
        -1
        >>> x = Variable(0)
        >>> -(x+1)
        -x0-1
        """
        return self*(-1)

    def __richcmp__(self, other, int op):
        """
        Construct :class:`Constraint`s

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> x+1 <  y-2
        -x0+x1-3>0
        >>> x+1 <= y-2
        -x0+x1-3>=0
        >>> x+1 == y-2
        x0-x1+3==0
        >>> x+1 >= y-2
        x0-x1+3>=0
        >>> x+1 >  y-2
        x0-x1+3>0
        """
        return _make_Constraint_from_richcmp(self, other, op)

    def __reduce__(self):
        """
        Pickle object

        Examples:

        >>> from ppl import Linear_Expression
        >>> le = loads(dumps(Linear_Expression([1,2,3],4)))
        >>> le.coefficients() == (1,2,3)
        True
        >>> le.inhomogeneous_term() == 4
        True
        """
        return (Linear_Expression, (self.coefficients(), self.inhomogeneous_term()))


cdef _wrap_Generator(PPL_Generator generator):
    """
    Wrap a C++ ``PPL_Generator`` into a Cython ``Generator``.
    """
    cdef Generator g = Generator(True)
    g.thisptr = new PPL_Generator(generator)
    return g



####################################################
cdef class Generator(object):
    r"""
    Wrapper for PPL's ``Generator`` class.

    An object of the class Generator is one of the following:

      * a line `\ell = (a_0, \dots, a_{n-1})^T`

      * a ray `r = (a_0, \dots, a_{n-1})^T`

      * a point `p = (\tfrac{a_0}{d}, \dots, \tfrac{a_{n-1}}{d})^T`

      * a closure point `c = (\tfrac{a_0}{d}, \dots, \tfrac{a_{n-1}}{d})^T`

    where `n` is the dimension of the space and, for points and
    closure points, `d` is the divisor.

    INPUT/OUTPUT:

    Use the helper functions :func:`line`, :func:`ray`, :func:`point`,
    and :func:`closure_point` to construct generators. Analogous class
    methods are also available, see :meth:`Generator.line`,
    :meth:`Generator.ray`, :meth:`Generator.point`,
    :meth:`Generator.closure_point`. Do not attempt to construct
    generators manually.

    .. NOTE::

        The generators are constructed from linear expressions. The
        inhomogeneous term is always silently discarded.

    Examples:

    >>> from ppl import Generator, Variable
    >>> x = Variable(0)
    >>> y = Variable(1)
    >>> Generator.line(5*x-2*y)
    line(5, -2)
    >>> Generator.ray(5*x-2*y)
    ray(5, -2)
    >>> Generator.point(5*x-2*y, 7)
    point(5/7, -2/7)
    >>> Generator.closure_point(5*x-2*y, 7)
    closure_point(5/7, -2/7)
    """
    def __cinit__(self, do_not_construct_manually=False):
        """
        The Cython constructor.

        See :class:`Variable` for documentation.

        Tests:

        >>> from ppl import Variable, line
        >>> x = Variable(0)
        >>> line(x)   # indirect doctest
        line(1)
        """
        assert(do_not_construct_manually)
        self.thisptr = NULL

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        assert self.thisptr!=NULL, 'Do not construct Generators manually!'
        del self.thisptr

    @classmethod
    def line(cls, expression):
        """
        Construct a line.

        INPUT:

        - ``expression`` -- a :class:`Linear_Expression` or something
          convertible to it (:class:`Variable` or integer).

        OUTPUT:

        A new :class:`Generator` representing the line.

        Raises a ``ValueError` if the homogeneous part of
        ``expression`` represents the origin of the vector space.

        Examples:

        >>> from ppl import Generator, Variable
        >>> y = Variable(1)
        >>> Generator.line(2*y)
        line(0, 1)
        >>> Generator.line(y)
        line(0, 1)
        >>> Generator.line(1)
        Traceback (most recent call last):
        ...
        ValueError: PPL::line(e):
        e == 0, but the origin cannot be a line.
        """
        cdef Linear_Expression e = Linear_Expression(expression)
        # This does not work as Cython gets confused by the private default ctor
        #   return _wrap_Generator(PPL_line(e.thisptr[0]))
        # workaround follows
        cdef Generator g = Generator(True)
        try:
            g.thisptr = new_line(e.thisptr[0])
        except BaseException:
            # g.thisptr must be set to something valid or g.__dealloc__() will segfault
            g.thisptr = new_point(e.thisptr[0],PPL_Coefficient(1))
            raise
        return g

    @classmethod
    def ray(cls, expression):
        """
        Construct a ray.

        INPUT:

        - ``expression`` -- a :class:`Linear_Expression` or something
          convertible to it (:class:`Variable` or integer).

        OUTPUT:

        A new :class:`Generator` representing the ray.

        Raises a ``ValueError` if the homogeneous part of
        ``expression`` represents the origin of the vector space.

        Examples:

        >>> from ppl import Generator, Variable
        >>> y = Variable(1)
        >>> Generator.ray(2*y)
        ray(0, 1)
        >>> Generator.ray(y)
        ray(0, 1)
        >>> Generator.ray(1)
        Traceback (most recent call last):
        ...
        ValueError: PPL::ray(e):
        e == 0, but the origin cannot be a ray.
        """
        cdef Linear_Expression e = Linear_Expression(expression)
        # This does not work as Cython gets confused by the private default ctor
        #   return _wrap_Generator(PPL_ray(e.thisptr[0]))
        # workaround follows
        cdef Generator g = Generator(True)
        try:
            g.thisptr = new_ray(e.thisptr[0])
        except BaseException:
            # g.thisptr must be set to something valid or g.__dealloc__() will segfault
            g.thisptr = new_point(e.thisptr[0],PPL_Coefficient(1))
            raise
        return g

    @classmethod
    def point(cls, expression=0, divisor=1):
        """
        Construct a point.

        INPUT:

        - ``expression`` -- a :class:`Linear_Expression` or something
          convertible to it (:class:`Variable` or integer).

        - ``divisor`` -- an integer.

        OUTPUT:

        A new :class:`Generator` representing the point.

        Raises a ``ValueError` if ``divisor==0``.

        Examples:

        >>> from ppl import Generator, Variable
        >>> y = Variable(1)
        >>> Generator.point(2*y+7, 3)
        point(0/3, 2/3)
        >>> Generator.point(y+7, 3)
        point(0/3, 1/3)
        >>> Generator.point(7, 3)
        point()
        >>> Generator.point(0, 0)
        Traceback (most recent call last):
        ...
        ValueError: PPL::point(e, d):
        d == 0.
        """
        cdef Linear_Expression e = Linear_Expression(expression)
        cdef long d = <long?> divisor
        # This does not work as Cython gets confused by the private default ctor
        #   return _wrap_Generator(PPL_point(e.thisptr[0], PPL_Coefficient(d.value)))
        # workaround follows
        cdef Generator g = Generator(True)
        try:
            g.thisptr = new_point(e.thisptr[0], PPL_Coefficient(d))
        except BaseException:
            # g.thisptr must be set to something valid or g.__dealloc__() will segfault
            g.thisptr = new_point(e.thisptr[0],PPL_Coefficient(1))
            raise
        return g

    @classmethod
    def closure_point(cls, expression=0, divisor=1):
        """
        Construct a closure point.

        A closure point is a point of the topological closure of a
        polyhedron that is not a point of the polyhedron itself.

        INPUT:

        - ``expression`` -- a :class:`Linear_Expression` or something
          convertible to it (:class:`Variable` or integer).

        - ``divisor`` -- an integer.

        OUTPUT:

        A new :class:`Generator` representing the point.

        Raises a ``ValueError` if ``divisor==0``.

        Examples:

        >>> from ppl import Generator, Variable
        >>> y = Variable(1)
        >>> Generator.closure_point(2*y+7, 3)
        closure_point(0/3, 2/3)
        >>> Generator.closure_point(y+7, 3)
        closure_point(0/3, 1/3)
        >>> Generator.closure_point(7, 3)
        closure_point()
        >>> Generator.closure_point(0, 0)
        Traceback (most recent call last):
        ...
        ValueError: PPL::closure_point(e, d):
        d == 0.
        """
        cdef Linear_Expression e = Linear_Expression(expression)
        cdef long d = <long> divisor
        # This does not work as Cython gets confused by the private default ctor
        #   return _wrap_Generator(PPL_closure_point(e.thisptr[0], PPL_Coefficient(d.value)))
        # workaround follows
        cdef Generator g = Generator(True)
        try:
            g.thisptr = new_closure_point(e.thisptr[0], PPL_Coefficient(d))
        except BaseException:
            # g.thisptr must be set to something valid or g.__dealloc__() will segfault
            g.thisptr = new_point(e.thisptr[0], PPL_Coefficient(1))
            raise
        return g

    def __repr__(self):
        """
        Return a string representation of the generator.

        OUTPUT:

        String.

        Examples:

        >>> from ppl import Generator, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 2*x-y+5
        >>> Generator.line(e)
        line(2, -1)
        >>> Generator.ray(e)
        ray(2, -1)
        >>> Generator.point(e, 3)
        point(2/3, -1/3)
        >>> Generator.closure_point(e, 3)
        closure_point(2/3, -1/3)
        """
        cdef PPL_GeneratorType t = self.thisptr.type()
        if t == LINE or t == RAY:
            div = ''
        elif t == POINT or t == CLOSURE_POINT:
            div = '/' + str(self.divisor())
        else:
            raise RuntimeError

        s = ', '.join(str(self.coefficient(Variable(i))) + div for i in range(self.space_dimension()))
        return PPL_GeneratorType_str(t) + '(' + s + ')'

    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUTPUT:

        Integer.

        Examples:

        >>> from ppl import Variable, point
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> point(x).space_dimension()
        1
        >>> point(y).space_dimension()
        2
        """
        return self.thisptr.space_dimension()

    def type(self):
        r"""
        Return the generator type of ``self``.

        OUTPUT:

        String. One of ``'line'``, ``'ray'``, ``'point'``, or
        ``'closure_point'``.

        Examples:

        >>> from ppl import Variable, point, closure_point, ray, line
        >>> x = Variable(0)
        >>> line(x).type()
        'line'
        >>> ray(x).type()
        'ray'
        >>> point(x,2).type()
        'point'
        >>> closure_point(x,2).type()
        'closure_point'
        """
        return PPL_GeneratorType_str(self.thisptr.type())

    def is_line(self):
        r"""
        Test whether ``self`` is a line.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import Variable, point, closure_point, ray, line
        >>> x = Variable(0)
        >>> line(x).is_line()
        True
        >>> ray(x).is_line()
        False
        >>> point(x,2).is_line()
        False
        >>> closure_point(x,2).is_line()
        False
        """
        return self.thisptr.is_line()

    def is_ray(self):
        r"""
        Test whether ``self`` is a ray.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import Variable, point, closure_point, ray, line
        >>> x = Variable(0)
        >>> line(x).is_ray()
        False
        >>> ray(x).is_ray()
        True
        >>> point(x,2).is_ray()
        False
        >>> closure_point(x,2).is_ray()
        False
        """
        return self.thisptr.is_ray()

    def is_line_or_ray(self):
        r"""
        Test whether ``self`` is a line or a ray.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import Variable, point, closure_point, ray, line
        >>> x = Variable(0)
        >>> line(x).is_line_or_ray()
        True
        >>> ray(x).is_line_or_ray()
        True
        >>> point(x,2).is_line_or_ray()
        False
        >>> closure_point(x,2).is_line_or_ray()
        False
        """
        return self.thisptr.is_line_or_ray()

    def is_point(self):
        r"""
        Test whether ``self`` is a point.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import Variable, point, closure_point, ray, line
        >>> x = Variable(0)
        >>> line(x).is_point()
        False
        >>> ray(x).is_point()
        False
        >>> point(x,2).is_point()
        True
        >>> closure_point(x,2).is_point()
        False
        """
        return self.thisptr.is_point()

    def is_closure_point(self):
        r"""
        Test whether ``self`` is a closure point.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import Variable, point, closure_point, ray, line
        >>> x = Variable(0)
        >>> line(x).is_closure_point()
        False
        >>> ray(x).is_closure_point()
        False
        >>> point(x,2).is_closure_point()
        False
        >>> closure_point(x,2).is_closure_point()
        True
        """
        return self.thisptr.is_closure_point()

    def coefficient(self, Variable v):
        """
        Return the coefficient of the variable ``v``.

        INPUT:

        - ``v`` -- a :class:`Variable`.

        OUTPUT:

        An integer.

        Examples:

        >>> from ppl import Variable, line
        >>> x = Variable(0)
        >>> line = line(3*x+1)
        >>> line
        line(1)
        >>> line.coefficient(x)
        1
        """
        return mpz_get_pyintlong(self.thisptr.coefficient(v.thisptr[0]).get_mpz_t())

    def coefficients(self):
        """
        Return the coefficients of the generator.

        See also :meth:`coefficient`.

        OUTPUT:

        A tuple of integers of length :meth:`space_dimension`.

        Examples:

            >>> from ppl import Variable, point
            >>> x = Variable(0);  y = Variable(1)
            >>> p = point(3*x+5*y+1, 2); p
            point(3/2, 5/2)
            >>> p.coefficients()
            (3, 5)
        """
        cdef int d = self.space_dimension()
        cdef int i
        coeffs = []
        for i in range(0,d):
            coeffs.append(mpz_get_pyintlong(self.thisptr.coefficient(PPL_Variable(i)).get_mpz_t()))
        return tuple(coeffs)

    def divisor(self):
        """
        If ``self`` is either a point or a closure point, return its
        divisor.

        OUTPUT:

        An integer. If ``self`` is a ray or a line, raises
        ``ValueError``.

        Examples:

        >>> from ppl import Generator, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> point = Generator.point(2*x-y+5)
        >>> point.divisor()
        1
        >>> line = Generator.line(2*x-y+5)
        >>> line.divisor()
        Traceback (most recent call last):
        ...
        ValueError: PPL::Generator::divisor():
        *this is neither a point nor a closure point.
        """
        return mpz_get_pyintlong(self.thisptr.divisor().get_mpz_t())

    def is_equivalent_to(self, Generator g):
        r"""
        Test whether ``self`` and ``g`` are equivalent.

        INPUT:

        - ``g`` -- a :class:`Generator`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` and ``g``
        are equivalent generators.

        Note that generators having different space dimensions are not
        equivalent.

        Examples:

        >>> from ppl import Generator, Variable, point, line
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> point(2*x    , 2).is_equivalent_to( point(x) )
        True
        >>> point(2*x+0*y, 2).is_equivalent_to( point(x) )
        False
        >>> line(4*x).is_equivalent_to(line(x))
        True
        """
        return self.thisptr.is_equivalent_to(g.thisptr[0])

    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        Examples:

        >>> cmd  = 'from ppl import Linear_Expression, Variable, point\n'
        >>> cmd += 'x = Variable(0)\n'
        >>> cmd += 'y = Variable(1)\n'
        >>> cmd += 'p = point(3*x+2*y)\n'
        >>> cmd += 'p.ascii_dump()\n'
        >>> import subprocess
        >>> proc = subprocess.Popen(['python', '-c', cmd], stderr=subprocess.PIPE)
        >>> out, err = proc.communicate()
        >>> print err
        size 3 1 3 2 P (C)
        <BLANKLINE>
        """
        self.thisptr.ascii_dump()

    def OK(self):
        """
        Check if all the invariants are satisfied.

        Examples:

        >>> from ppl import Linear_Expression, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> e = 3*x+2*y+1
        >>> e.OK()
        True
    """
        return self.thisptr.OK()

    def __reduce__(self):
        """
        Pickle object.

        Tests:

        >>> from ppl import Generator, Variable, line, ray, point, closure_point
        >>> from pickle import loads,dumps
        >>> x = Variable(0); y = Variable(1)
        >>> loads(dumps(Generator.point(2*x+7*y, 3)))
        point(2/3, 7/3)
        >>> loads(dumps(Generator.closure_point(2*x+7*y, 3)))
        closure_point(2/3, 7/3)
        >>> loads(dumps(Generator.line(2*x+7*y)))
        line(2, 7)
        >>> loads(dumps(Generator.ray(2*x+7*y)))
        ray(2, 7)
        """
        t = self.thisptr.type()
        le = Linear_Expression(self.coefficients(), 0)
        if t == LINE:
            return (line, (le,))
        elif t == RAY:
            return (ray, (le,))
        elif t == POINT:
            return (point, (le, self.divisor()))
        elif t == CLOSURE_POINT:
            return (closure_point, (le, self.divisor()))
        else:
            raise RuntimeError


#############
# shortcuts #
#############
line          = Generator.line
ray           = Generator.ray
point         = Generator.point
closure_point = Generator.closure_point


####################################################
### Generator_System  ##############################
####################################################
####################################################
cdef class Generator_System(_mutable_or_immutable):
    """
    Wrapper for PPL's ``Generator_System`` class.

    An object of the class Generator_System is a system of generators,
    i.e., a multiset of objects of the class Generator (lines, rays,
    points and closure points). When inserting generators in a system,
    space dimensions are automatically adjusted so that all the
    generators in the system are defined on the same vector space. A
    system of generators which is meant to define a non-empty
    polyhedron must include at least one point: the reason is that
    lines, rays and closure points need a supporting point (lines and
    rays only specify directions while closure points only specify
    points in the topological closure of the NNC polyhedron).

    Examples:

        >>> from ppl import Generator_System, Variable, line, ray, point, closure_point
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> gs = Generator_System(line(5*x-2*y))
        >>> gs.insert(ray(6*x-3*y))
        >>> gs.insert(point(2*x-7*y, 5))
        >>> gs.insert(closure_point(9*x-1*y, 2))
        >>> gs
        Generator_System {line(5, -2), ray(2, -1), point(2/5, -7/5), closure_point(9/2, -1/2)}
    """
    def __cinit__(self, arg=None):
        """
        The Cython constructor.

        See :class:`Generator_System` for documentation.

        Tests:

        >>> from ppl import Generator_System
        >>> Generator_System()   # indirect doctest
        Generator_System {}
        """
        if arg is None:
            self.thisptr = new PPL_Generator_System()
            return
        if isinstance(arg, Generator):
            g = <Generator>arg
            self.thisptr = new PPL_Generator_System(g.thisptr[0])
            return
        if isinstance(arg, Generator_System):
            gs = <Generator_System>arg
            self.thisptr = new PPL_Generator_System(gs.thisptr[0])
            return
        if isinstance(arg, (list,tuple)):
            self.thisptr = new PPL_Generator_System()
            for generator in arg:
                self.insert(generator)
            return
        raise ValueError, 'Cannot initialize with '+str(arg)+'.'

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr

    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUTPUT:

        Integer.

        Examples:

        >>> from ppl import Variable, Generator_System, point
        >>> x = Variable(0)
        >>> gs = Generator_System( point(3*x) )
        >>> gs.space_dimension()
        1
        """
        return self.thisptr.space_dimension()

    def clear(self):
        r"""
        Removes all generators from the generator system and sets its
        space dimension to 0.

        Examples:

        >>> from ppl import Variable, Generator_System, point
        >>> x = Variable(0)
        >>> gs = Generator_System( point(3*x) ); gs
        Generator_System {point(3/1)}
        >>> gs.clear()
        >>> gs
        Generator_System {}
        """
        self.assert_mutable('The Generator_System is not mutable!')
        self.thisptr.clear()

    def insert(self, Generator g):
        """
        Insert ``g`` into the generator system.

        The number of space dimensions of ``self`` is increased, if needed.

        INPUT:

        - ``g`` -- a :class:`Generator`.

        Examples:

        >>> from ppl import Variable, Generator_System, point
        >>> x = Variable(0)
        >>> gs = Generator_System( point(3*x) )
        >>> gs.insert( point(-3*x) )
        >>> gs
        Generator_System {point(3/1), point(-3/1)}
        """
        self.assert_mutable('The Generator_System is not mutable!')
        self.thisptr.insert(g.thisptr[0])

    def empty(self):
        """
        Return ``True`` if and only if ``self`` has no generators.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import Variable, Generator_System, point
        >>> x = Variable(0)
        >>> gs = Generator_System()
        >>> gs.empty()
        True
        >>> gs.insert( point(-3*x) )
        >>> gs.empty()
        False
        """
        return self.thisptr.empty()

    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        Examples:

        >>> cmd  = 'from ppl import Generator_System, point, Variable\n'
        >>> cmd += 'x = Variable(0)\n'
        >>> cmd += 'y = Variable(1)\n'
        >>> cmd += 'gs = Generator_System( point(3*x+2*y+1) )\n'
        >>> cmd += 'gs.ascii_dump()\n'
        >>> import subprocess
        >>> proc = subprocess.Popen(['python', '-c', cmd], stderr=subprocess.PIPE)
        >>> out, err = proc.communicate()
        >>> print err
        topology NECESSARILY_CLOSED
        1 x 2 SPARSE (sorted)
        index_first_pending 1
        size 3 1 3 2 P (C)
        <BLANKLINE>
        """
        self.thisptr.ascii_dump()

    def OK(self):
        """
        Check if all the invariants are satisfied.

        Examples:

        >>> from ppl import Variable, Generator_System, point
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> gs = Generator_System( point(3*x+2*y+1) )
        >>> gs.OK()
        True
        """
        return self.thisptr.OK()

    def __len__(self):
        """
        Return the number of generators in the system.

        >>> from ppl import Variable, Generator_System, point
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> gs = Generator_System()
        >>> gs.insert(point(3*x+2*y))
        >>> gs.insert(point(x))
        >>> gs.insert(point(y))
        >>> len(gs)
        3
        """
        return sum([1 for g in self])

    def __iter__(self):
        """
        Iterate through the generators of the system.

        Examples:

        >>> from ppl import Generator_System, Variable, point
        >>> x = Variable(0)
        >>> gs = Generator_System(point(3*x))
        >>> iter = gs.__iter__()
        >>> next(iter)
        point(3/1)
        """
        return Generator_System_iterator(self)

    def __getitem__(self, int k):
        """
        Return the ``k``-th generator.

        The correct way to read the individual generators is to
        iterate over the generator system. This method is for
        convenience only.

        INPUT:

        - ``k`` -- integer. The index of the generator.

        OUTPUT:

        The ``k``-th constraint of the generator system.

        Examples:

        >>> from ppl import Generator_System, Variable, point
        >>> x = Variable(0)
        >>> gs = Generator_System()
        >>> gs.insert(point(3*x))
        >>> gs.insert(point(-2*x))
        >>> gs
        Generator_System {point(3/1), point(-2/1)}
        >>> gs[0]
        point(3/1)
        >>> gs[1]
        point(-2/1)
        """
        if k < 0:
            raise IndexError('index must be nonnegative')
        iterator = iter(self)
        try:
            for i in range(k):
                next(iterator)
        except StopIteration:
            raise IndexError('index is past-the-end')
        return next(iterator)

    def __repr__(self):
        r"""
        Return a string representation of the generator system.

        OUTPUT:

        A string.

        Examples:

        >>> from ppl import Generator_System, Variable, point, ray
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> gs = Generator_System(point(3*x+2*y+1))
        >>> gs.insert(ray(x))
        >>> gs.__repr__()
        'Generator_System {point(3/1, 2/1), ray(1, 0)}'
        """
        s = 'Generator_System {'
        s += ', '.join([ repr(g) for g in self ])
        s += '}'
        return s

    def __reduce__(self):
        """
        Pickle object.

        Tests:

        >>> from ppl import Generator_System, Variable, point, ray
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> gs = Generator_System((point(3*x+2*y+1), ray(x)));  gs
        Generator_System {point(3/1, 2/1), ray(1, 0)}
        >>> loads(dumps(gs))
        Generator_System {point(3/1, 2/1), ray(1, 0)}
        """
        return (Generator_System, (tuple(self), ))



####################################################
### Generator_System_iterator ######################
####################################################

####################################################
cdef class Generator_System_iterator(object):
    """
    Wrapper for PPL's ``Generator_System::const_iterator`` class.

    Examples:

    >>> from ppl import Generator_System, Variable, line, ray, point, closure_point, Generator_System_iterator
    >>> x = Variable(0)
    >>> y = Variable(1)
    >>> gs = Generator_System( line(5*x-2*y) )
    >>> gs.insert( ray(6*x-3*y) )
    >>> gs.insert( point(2*x-7*y, 5) )
    >>> gs.insert( closure_point(9*x-1*y, 2) )
    >>> next(Generator_System_iterator(gs))
    line(5, -2)
    >>> list(gs)
    [line(5, -2), ray(2, -1), point(2/5, -7/5), closure_point(9/2, -1/2)]
    """
    def __cinit__(self, Generator_System gs):
        r"""
        The Cython constructor.

        Tests:

        >>> from ppl import Generator_System, Generator_System_iterator
        >>> iter = Generator_System_iterator(Generator_System())   # indirect doctest
        """
        self.gs = gs
        self.gsi_ptr = init_gs_iterator(gs.thisptr[0])

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        delete_gs_iterator(self.gsi_ptr)

    def __next__(Generator_System_iterator self):
        r"""
        The next iteration.

        OUTPUT:

        A :class:`Generator`.

        Examples:

        >>> from ppl import Generator_System, Variable, point, Generator_System_iterator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> gs = Generator_System( point(5*x-2*y) )
        >>> next(Generator_System_iterator(gs))
        point(5/1, -2/1)
        """
        if is_end_gs_iterator((<Generator_System>self.gs).thisptr[0], self.gsi_ptr):
            raise StopIteration
        return _wrap_Generator(next_gs_iterator(self.gsi_ptr))


####################################################
### Constraint ######################################
####################################################

####################################################
cdef _make_Constraint_from_richcmp(lhs_, rhs_, op):
    cdef Linear_Expression lhs = Linear_Expression(lhs_)
    cdef Linear_Expression rhs = Linear_Expression(rhs_)
    if op==0:      # <   0
        return _wrap_Constraint(lhs.thisptr[0] <  rhs.thisptr[0])
    elif op==1:    # <=  1
        return _wrap_Constraint(lhs.thisptr[0] <= rhs.thisptr[0])
    elif op==2:    # ==  2
        return _wrap_Constraint(lhs.thisptr[0] == rhs.thisptr[0])
    elif op==4:    # >   4
        return _wrap_Constraint(lhs.thisptr[0] >  rhs.thisptr[0])
    elif op==5:    # >=  5
        return _wrap_Constraint(lhs.thisptr[0] >= rhs.thisptr[0])
    elif op==3:    # !=  3
        raise NotImplementedError
    else:
        assert(False)


####################################################
cdef class Constraint(object):
    """
    Wrapper for PPL's ``Constraint`` class.

    An object of the class ``Constraint`` is either:

    * an equality `\sum_{i=0}^{n-1} a_i x_i + b = 0`

    * a non-strict inequality `\sum_{i=0}^{n-1} a_i x_i + b \geq 0`

    * a strict inequality `\sum_{i=0}^{n-1} a_i x_i + b > 0`

    where `n` is the dimension of the space, `a_i` is the integer
    coefficient of variable `x_i`, and `b_i` is the integer
    inhomogeneous term.

    INPUT/OUTPUT:

    You construct constraints by writing inequalities in
    :class:`Linear_Expression`. Do not attempt to manually construct
    constraints.

    Examples:

    >>> from ppl import Constraint, Variable, Linear_Expression
    >>> x = Variable(0)
    >>> y = Variable(1)
    >>> 5*x-2*y >  x+y-1
    4*x0-3*x1+1>0
    >>> 5*x-2*y >= x+y-1
    4*x0-3*x1+1>=0
    >>> 5*x-2*y == x+y-1
    4*x0-3*x1+1==0
    >>> 5*x-2*y <= x+y-1
    -4*x0+3*x1-1>=0
    >>> 5*x-2*y <  x+y-1
    -4*x0+3*x1-1>0
    >>> x > 0
    x0>0

    Special care is needed if the left hand side is a constant:

    >>> 0 == 1    # watch out!
    False
    >>> Linear_Expression(0) == 1
    -1==0
    """
    def __cinit__(self, do_not_construct_manually=False):
        """
        The Cython constructor.

        See :class:`Constraint` for documentation.

        Tests:

            >>> from ppl import Constraint, Variable, Linear_Expression
            >>> x = Variable(0)
            >>> x>0   # indirect doctest
            x0>0
        """
        assert(do_not_construct_manually)
        self.thisptr = NULL

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        assert self.thisptr!=NULL, 'Do not construct Constraints manually!'
        del self.thisptr

    def __repr__(self):
        """
        Return a string representation of the constraint.

        OUTPUT:

        String.

        Examples:

        >>> from ppl import Constraint, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> (2*x-y+5 >  x).__repr__()
        'x0-x1+5>0'
        >>> (2*x-y+5 == x).__repr__()
        'x0-x1+5==0'
        >>> (2*x-y+5 >= x).__repr__()
        'x0-x1+5>=0'
        """
        e = sum(self.coefficient(x)*x
                for x in (Variable(i)
                          for i in range(self.space_dimension())))
        e += self.inhomogeneous_term()
        s = repr(e)
        t = self.thisptr.type()
        if t == EQUALITY:
            s += '==0'
        elif t == NONSTRICT_INEQUALITY:
            s += '>=0'
        elif t == STRICT_INEQUALITY:
            s += '>0'
        else:
            raise RuntimeError
        return s

    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUTPUT:

        Integer.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> (x>=0).space_dimension()
        1
        >>> (y==1).space_dimension()
        2
        """
        return self.thisptr.space_dimension()

    def type(self):
        r"""
        Return the constraint type of ``self``.

        OUTPUT:

        String. One of ``'equality'``, ``'nonstrict_inequality'``, or
        ``'strict_inequality'``.

        Examples:

            >>> from ppl import Variable
            >>> x = Variable(0)
            >>> (x==0).type()
            'equality'
            >>> (x>=0).type()
            'nonstrict_inequality'
            >>> (x>0).type()
            'strict_inequality'
        """
        t = self.thisptr.type()
        if t == EQUALITY:
            return 'equality'
        elif t == NONSTRICT_INEQUALITY:
            return 'nonstrict_inequality'
        elif t == STRICT_INEQUALITY:
            return 'strict_inequality'
        else:
            raise RuntimeError

    def is_equality(self):
        r"""
        Test whether ``self`` is an equality.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is an
        equality constraint.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> (x==0).is_equality()
        True
        >>> (x>=0).is_equality()
        False
        >>> (x>0).is_equality()
        False
        """
        return self.thisptr.is_equality()

    def is_inequality(self):
        r"""
        Test whether ``self`` is an inequality.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is an
        inequality constraint, either strict or non-strict.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> (x==0).is_inequality()
        False
        >>> (x>=0).is_inequality()
        True
        >>> (x>0).is_inequality()
        True
        """
        return self.thisptr.is_inequality()

    def is_nonstrict_inequality(self):
        r"""
        Test whether ``self`` is a non-strict inequality.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is an
        non-strict inequality constraint.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> (x==0).is_nonstrict_inequality()
        False
        >>> (x>=0).is_nonstrict_inequality()
        True
        >>> (x>0).is_nonstrict_inequality()
        False
        """
        return self.thisptr.is_nonstrict_inequality()

    def is_strict_inequality(self):
        r"""
        Test whether ``self`` is a strict inequality.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is an
        strict inequality constraint.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> (x==0).is_strict_inequality()
        False
        >>> (x>=0).is_strict_inequality()
        False
        >>> (x>0).is_strict_inequality()
        True
        """
        return self.thisptr.is_strict_inequality()

    def coefficient(self, Variable v):
        """
        Return the coefficient of the variable ``v``.

        INPUT:

        - ``v`` -- a :class:`Variable`.

        OUTPUT:

        An integer.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> ineq = 3*x+1 > 0
        >>> ineq.coefficient(x)
        3
        >>> y = Variable(1)
        >>> ineq = 3**50 * y + 2 > 1
        >>> ineq.coefficient(y)
        717897987691852588770249L
        >>> ineq.coefficient(x)
        0
        """
        return mpz_get_pyintlong(self.thisptr.coefficient(v.thisptr[0]).get_mpz_t())

    def coefficients(self):
        """
        Return the coefficients of the constraint.

        See also :meth:`coefficient`.

        OUTPUT:

        A tuple of integers of length :meth:`space_dimension`.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0);  y = Variable(1)
        >>> ineq = ( 3*x+5*y+1 ==  2);  ineq
        3*x0+5*x1-1==0
        >>> ineq.coefficients()
        (3, 5)
        """
        cdef int d = self.space_dimension()
        cdef int i
        coeffs = []
        for i in range(0,d):
            coeffs.append(mpz_get_pyintlong(self.thisptr.coefficient(PPL_Variable(i)).get_mpz_t()))
        return tuple(coeffs)

    def inhomogeneous_term(self):
        """
        Return the inhomogeneous term of the constraint.

        OUTPUT:

        Integer.

        Examples:

        >>> from ppl import Variable
        >>> y = Variable(1)
        >>> ineq = 10+y > 9
        >>> ineq
        x1+1>0
        >>> ineq.inhomogeneous_term()
        1
        >>> ineq = 2**66 + y > 0
        >>> ineq.inhomogeneous_term()
        73786976294838206464L
        """
        return mpz_get_pyintlong(self.thisptr.inhomogeneous_term().get_mpz_t())

    def is_tautological(self):
        r"""
        Test whether ``self`` is a tautological constraint.

        A tautology can have either one of the following forms:

        * an equality: `\sum 0 x_i + 0 = 0`,

        * a non-strict inequality: `\sum 0 x_i + b \geq 0` with `b\geq 0`, or

        * a strict inequality: `\sum 0 x_i + b > 0` with `b> 0`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is a
        tautological constraint.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> (x==0).is_tautological()
        False
        >>> (0*x>=0).is_tautological()
        True
        """
        return self.thisptr.is_tautological()

    def is_inconsistent(self):
        r"""
        Test whether ``self`` is an inconsistent constraint, that is, always false.

        An inconsistent constraint can have either one of the
        following forms:

        * an equality: `\sum 0 x_i + b = 0` with `b\not=0`,

        * a non-strict inequality: `\sum 0 x_i + b \geq 0` with `b< 0`, or

        * a strict inequality: `\sum 0 x_i + b > 0` with `b\leq 0`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` is an
        inconsistent constraint.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> (x==1).is_inconsistent()
        False
        >>> (0*x>=1).is_inconsistent()
        True
        """
        return self.thisptr.is_inconsistent()

    def is_equivalent_to(self, Constraint c):
        r"""
        Test whether ``self`` and ``c`` are equivalent.

        INPUT:

        - ``c`` -- a :class:`Constraint`.

        OUTPUT:

        Boolean. Returns ``True`` if and only if ``self`` and ``c``
        are equivalent constraints.

        Note that constraints having different space dimensions are
        not equivalent. However, constraints having different types
        may nonetheless be equivalent, if they both are tautologies or
        inconsistent.

        Examples:

        >>> from ppl import Variable, Linear_Expression
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> (x > 0).is_equivalent_to(Linear_Expression(0) < x)
        True
        >>> (x > 0).is_equivalent_to(0*y < x)
        False
        >>> (0*x > 1).is_equivalent_to(0*x == -2)
        True
        """
        return self.thisptr.is_equivalent_to(c.thisptr[0])

    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        Examples:

        >>> cmd  = 'from ppl import Linear_Expression, Variable\n'
        >>> cmd += 'x = Variable(0)\n'
        >>> cmd += 'y = Variable(1)\n'
        >>> cmd += 'e = (3*x+2*y+1 > 0)\n'
        >>> cmd += 'e.ascii_dump()\n'
        >>> import subprocess
        >>> proc = subprocess.Popen(['python', '-c', cmd], stderr=subprocess.PIPE)
        >>> out, err = proc.communicate()
        >>> print err
        size 4 1 3 2 -1 > (NNC)
        <BLANKLINE>
        """
        self.thisptr.ascii_dump()

    def OK(self):
        """
        Check if all the invariants are satisfied.

        Examples:

        >>> from ppl import Linear_Expression, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> ineq = (3*x+2*y+1>=0)
        >>> ineq.OK()
        True
        """
        return self.thisptr.OK()

    def __reduce__(self):
        """
        Pickle object.

        Examples:

        >>> from ppl import Linear_Expression, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> loads(dumps(3*x+2*y+1>=5))
        3*x0+2*x1-4>=0
        >>> loads(dumps(3*x+2*y+1>5))
        3*x0+2*x1-4>0
        >>> loads(dumps(3*x+2*y+1==5))
        3*x0+2*x1-4==0
        """
        le = Linear_Expression(self.coefficients(), self.inhomogeneous_term())
        if self.is_nonstrict_inequality():
            return (inequality, (le, ))
        elif self.is_strict_inequality():
            return (strict_inequality, (le, ))
        elif self.is_equality():
            return (equation, (le, ))
        else:
            raise RuntimeError


####################################################
def inequality(expression):
    """
    Constuct an inequality.

    INPUT:

    - ``expression`` -- a :class:`Linear_Expression`.

    OUTPUT:

    The inequality ``expression`` >= 0.

    Examples:

    >>> from ppl import Variable, inequality
    >>> y = Variable(1)
    >>> 2*y+1 >= 0
    2*x1+1>=0
    >>> inequality(2*y+1)
    2*x1+1>=0
    """
    return expression >= 0


####################################################
def strict_inequality(expression):
    """
    Constuct a strict inequality.

    INPUT:

    - ``expression`` -- a :class:`Linear_Expression`.

    OUTPUT:

    The inequality ``expression`` > 0.

    Examples:

    >>> from ppl import Variable, strict_inequality
    >>> y = Variable(1)
    >>> 2*y+1 > 0
    2*x1+1>0
    >>> strict_inequality(2*y+1)
    2*x1+1>0
    """
    return expression > 0


####################################################
def equation(expression):
    """
    Constuct an equation.

    INPUT:

    - ``expression`` -- a :class:`Linear_Expression`.

    OUTPUT:

    The equation ``expression`` == 0.

    Examples:

    >>> from ppl import Variable, equation
    >>> y = Variable(1)
    >>> 2*y+1 == 0
    2*x1+1==0
    >>> equation(2*y+1)
    2*x1+1==0
    """
    return expression == 0



####################################################
### Constraint_System  ##############################
####################################################


####################################################
cdef class Constraint_System(object):
    """
    Wrapper for PPL's ``Constraint_System`` class.

    An object of the class Constraint_System is a system of
    constraints, i.e., a multiset of objects of the class
    Constraint. When inserting constraints in a system, space
    dimensions are automatically adjusted so that all the constraints
    in the system are defined on the same vector space.

    Examples:

        >>> from ppl import Constraint_System, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System( 5*x-2*y > 0 )
        >>> cs.insert( 6*x<3*y )
        >>> cs.insert( x >= 2*x-7*y )
        >>> cs
        Constraint_System {5*x0-2*x1>0, -2*x0+x1>0, -x0+7*x1>=0}
    """
    def __cinit__(self, arg=None):
        """
        The Cython constructor.

        See :class:`Constraint_System` for documentation.

        Tests:

        >>> from ppl import Constraint_System
        >>> Constraint_System()
        Constraint_System {}
        """
        if arg is None:
            self.thisptr = new PPL_Constraint_System()
        elif isinstance(arg, Constraint):
            g = <Constraint>arg
            self.thisptr = new PPL_Constraint_System(g.thisptr[0])
        elif isinstance(arg, Constraint_System):
            gs = <Constraint_System>arg
            self.thisptr = new PPL_Constraint_System(gs.thisptr[0])
        elif isinstance(arg, (list,tuple)):
            self.thisptr = new PPL_Constraint_System()
            for constraint in arg:
                self.insert(constraint)
        else:
            raise TypeError('cannot initialize from {!r}'.format(arg))

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr

    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUTPUT:

        Integer.

        Examples:

        >>> from ppl import Variable, Constraint_System
        >>> x = Variable(0)
        >>> cs = Constraint_System( x>0 )
        >>> cs.space_dimension()
        1
        """
        return self.thisptr.space_dimension()

    def has_equalities(self):
        r"""
        Tests whether ``self`` contains one or more equality constraints.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import Variable, Constraint_System
        >>> x = Variable(0)
        >>> cs = Constraint_System()
        >>> cs.insert( x>0 )
        >>> cs.insert( x<0 )
        >>> cs.has_equalities()
        False
        >>> cs.insert( x==0 )
        >>> cs.has_equalities()
        True
        """
        return self.thisptr.has_equalities()

    def has_strict_inequalities(self):
        r"""
        Tests whether ``self`` contains one or more strict inequality constraints.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import Variable, Constraint_System
        >>> x = Variable(0)
        >>> cs = Constraint_System()
        >>> cs.insert( x>=0 )
        >>> cs.insert( x==-1 )
        >>> cs.has_strict_inequalities()
        False
        >>> cs.insert( x>0 )
        >>> cs.has_strict_inequalities()
        True
        """
        return self.thisptr.has_strict_inequalities()

    def clear(self):
        r"""
        Removes all constraints from the constraint system and sets its
        space dimension to 0.

        Examples:

        >>> from ppl import Variable, Constraint_System
        >>> x = Variable(0)
        >>> cs = Constraint_System(x>0)
        >>> cs
        Constraint_System {x0>0}
        >>> cs.clear()
        >>> cs
        Constraint_System {}
        """
        self.assert_mutable('The Constraint_System is not mutable!')
        self.thisptr.clear()

    def insert(self, Constraint c):
        """
        Insert ``c`` into the constraint system.

        INPUT:

        - ``c`` -- a :class:`Constraint`.

        Examples:

        >>> from ppl import Variable, Constraint_System
        >>> x = Variable(0)
        >>> cs = Constraint_System()
        >>> cs.insert( x>0 )
        >>> cs
        Constraint_System {x0>0}
        """
        self.assert_mutable('The Constraint_System is not mutable!')
        self.thisptr.insert(c.thisptr[0])

    def empty(self):
        """
        Return ``True`` if and only if ``self`` has no constraints.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import Variable, Constraint_System, point
        >>> x = Variable(0)
        >>> cs = Constraint_System()
        >>> cs.empty()
        True
        >>> cs.insert( x>0 )
        >>> cs.empty()
        False
        """
        return self.thisptr.empty()

    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        Examples:

        >>> cmd  = 'from ppl import Constraint_System, Variable\n'
        >>> cmd += 'x = Variable(0)\n'
        >>> cmd += 'y = Variable(1)\n'
        >>> cmd += 'cs = Constraint_System( 3*x > 2*y+1 )\n'
        >>> cmd += 'cs.ascii_dump()\n'
        >>> import subprocess
        >>> proc = subprocess.Popen(['python', '-c', cmd], stderr=subprocess.PIPE)
        >>> out, err = proc.communicate()
        >>> print err
        topology NOT_NECESSARILY_CLOSED
        1 x 2 SPARSE (sorted)
        index_first_pending 1
        size 4 -1 3 -2 -1 > (NNC)
        <BLANKLINE>
        """
        self.thisptr.ascii_dump()

    def OK(self):
        """
        Check if all the invariants are satisfied.

        Examples:

        >>> from ppl import Variable, Constraint_System
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System( 3*x+2*y+1 <= 10 )
        >>> cs.OK()
        True
        """
        return self.thisptr.OK()

    def __len__(self):
        """
        Return the number of constraints in the system.

        Examples:

        >>> from ppl import Variable, Constraint_System
        >>> x = Variable(0)
        >>> cs = Constraint_System( x>0 )
        >>> cs.insert( x<1 )
        >>> len(cs)
        2
        """
        return sum(1 for c in self)

    def __iter__(self):
        """
        Iterate through the constraints of the system.

        Examples:

        >>> from ppl import Variable, Constraint_System
        >>> x = Variable(0)
        >>> cs = Constraint_System( x>0 )
        >>> iter = cs.__iter__()
        >>> next(iter)
        x0>0
        >>> list(cs)   # uses __iter__() internally
        [x0>0]
        """
        return Constraint_System_iterator(self)

    def __getitem__(self, int k):
        """
        Return the k-th constraint.

        The correct way to read the individual constraints is to
        iterate over the constraint system. This method is for
        convenience only.

        INPUT:

        - ``k`` -- integer. The index of the constraint.

        OUTPUT:

        The `k`-th constraint of the constraint system.

        Examples:

        >>> from ppl import Variable, Constraint_System
        >>> x = Variable(0)
        >>> cs = Constraint_System( x>0 )
        >>> cs.insert( x<1 )
        >>> cs
        Constraint_System {x0>0, -x0+1>0}
        >>> cs[0]
        x0>0
        >>> cs[1]
        -x0+1>0
        """
        if k < 0:
            raise IndexError('index must be nonnegative')
        iterator = iter(self)
        try:
            for i in range(k):
                next(iterator)
        except StopIteration:
            raise IndexError('index is past-the-end')
        return next(iterator)

    def __repr__(self):
        r"""
        Return a string representation of the constraint system.

        OUTPUT:

        A string.

        Examples:

        >>> from ppl import Constraint_System, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System([3*x+2*y+1 < 3, 0*x>x+1])
        >>> cs.__repr__()
        'Constraint_System {-3*x0-2*x1+2>0, -x0-1>0}'
        """
        s = 'Constraint_System {'
        s += ', '.join([ repr(c) for c in self ])
        s += '}'
        return s

    def __reduce__(self):
        """
        Pickle object.

        Tests:

        >>> from ppl import Constraint_System, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System([3*x+2*y+1 < 3, 0*x>x+1]);  cs
        Constraint_System {-3*x0-2*x1+2>0, -x0-1>0}
        >>> loads(dumps(cs))
        Constraint_System {-3*x0-2*x1+2>0, -x0-1>0}
        """
        return (Constraint_System, (tuple(self), ))


####################################################
### Constraint_System_iterator #####################
####################################################

####################################################
cdef class Constraint_System_iterator(object):
    """
    Wrapper for PPL's ``Constraint_System::const_iterator`` class.

    Examples:

        >>> from ppl import Constraint_System, Variable, Constraint_System_iterator
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System( 5*x < 2*y )
        >>> cs.insert( 6*x-3*y==0 )
        >>> cs.insert( x >= 2*x-7*y )
        >>> next(Constraint_System_iterator(cs))
        -5*x0+2*x1>0
        >>> list(cs)
        [-5*x0+2*x1>0, 2*x0-x1==0, -x0+7*x1>=0]
    """
    def __cinit__(self, Constraint_System cs):
        """
        The Cython constructor.

        See :class:`Constraint_System_iterator` for documentation.

        Tests:

        >>> from ppl import Constraint_System, Constraint_System_iterator
        >>> iter = Constraint_System_iterator( Constraint_System() )   # indirect doctest
        """
        self.cs = cs
        self.csi_ptr = init_cs_iterator(cs.thisptr[0])

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        delete_cs_iterator(self.csi_ptr)

    def __next__(Constraint_System_iterator self):
        r"""
        The next iteration.

        OUTPUT:

        A :class:`Generator`.

        Examples:

        >>> from ppl import Constraint_System, Variable, Constraint_System_iterator
        >>> x = Variable(0)
        >>> cs = Constraint_System( 5*x > 0 )
        >>> next(Constraint_System_iterator(cs))
        x0>0
        """
        if is_end_cs_iterator((<Constraint_System>self.cs).thisptr[0], self.csi_ptr):
            raise StopIteration
        return _wrap_Constraint(next_cs_iterator(self.csi_ptr))


cdef class Poly_Gen_Relation(object):
    r"""
    Wrapper for PPL's ``Poly_Con_Relation`` class.

    INPUT/OUTPUT:

    You must not construct :class:`Poly_Gen_Relation` objects
    manually. You will usually get them from
    :meth:`~sage.libs.ppl.Polyhedron.relation_with`. You can also get
    pre-defined relations from the class methods :meth:`nothing` and
    :meth:`subsumes`.

    Examples:

    >>> from ppl import Poly_Gen_Relation
    >>> nothing = Poly_Gen_Relation.nothing(); nothing
    nothing
    >>> subsumes = Poly_Gen_Relation.subsumes(); subsumes
    subsumes
    >>> nothing.implies( subsumes )
    False
    >>> subsumes.implies( nothing )
    True
    """
    def __cinit__(self, do_not_construct_manually=False):
        """
        The Cython constructor.

        See :class:`Poly_Gen_Relation` for documentation.

        Tests:

        >>> from ppl import Poly_Gen_Relation
        >>> Poly_Gen_Relation.nothing()
        nothing
        """
        assert(do_not_construct_manually)
        self.thisptr = NULL

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        assert self.thisptr!=NULL, 'Do not construct Poly_Gen_Relation objects manually!'
        del self.thisptr

    def implies(self, Poly_Gen_Relation y):
        r"""
        Test whether ``self`` implies ``y``.

        INPUT:

        - ``y`` -- a :class:`Poly_Gen_Relation`.

        OUTPUT:

        Boolean. ``True`` if and only if ``self`` implies ``y``.

        Examples:

        >>> from ppl import Poly_Gen_Relation
        >>> nothing = Poly_Gen_Relation.nothing()
        >>> nothing.implies( nothing )
        True
        """
        return self.thisptr.implies(y.thisptr[0])

    @classmethod
    def nothing(cls):
        r"""
        Return the assertion that says nothing.

        OUTPUT:

        A :class:`Poly_Gen_Relation`.

        Examples:

        >>> from ppl import Poly_Gen_Relation
        >>> Poly_Gen_Relation.nothing()
        nothing
        """
        return _wrap_Poly_Gen_Relation(PPL_Poly_Gen_Relation_nothing())

    @classmethod
    def subsumes(cls):
        r"""
        Return the assertion "Adding the generator would not change
        the polyhedron".

        OUTPUT:

        A :class:`Poly_Gen_Relation`.

        Examples:

        >>> from ppl import Poly_Gen_Relation
        >>> Poly_Gen_Relation.subsumes()
        subsumes
        """
        return _wrap_Poly_Gen_Relation(PPL_Poly_Gen_Relation_subsumes())

    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        Examples:

        >>> cmd  = 'from ppl import Poly_Gen_Relation\n'
        >>> cmd += 'Poly_Gen_Relation.nothing().ascii_dump()\n'
        >>> import subprocess
        >>> proc = subprocess.Popen(['python', '-c', cmd], stderr=subprocess.PIPE)
        >>> out, err = proc.communicate()
        >>> print err
        NOTHING
        >>> proc.returncode
        0
        """
        self.thisptr.ascii_dump()

    def OK(self, check_non_empty=False):
        """
        Check if all the invariants are satisfied.

        Examples:

            >>> from ppl import Poly_Gen_Relation
            >>> Poly_Gen_Relation.nothing().OK()
            True
        """
        return self.thisptr.OK()

    def __repr__(self):
        r"""
        Return a string representation.

        OUTPUT:

        String.

        Examples:

        >>> from ppl import Poly_Gen_Relation
        >>> Poly_Gen_Relation.nothing().__repr__()
        'nothing'
        """
        if self.implies(Poly_Gen_Relation.subsumes()):
            return 'subsumes'
        else:
            return 'nothing'




cdef class Poly_Con_Relation(object):
    r"""
    Wrapper for PPL's ``Poly_Con_Relation`` class.

    INPUT/OUTPUT:

    You must not construct :class:`Poly_Con_Relation` objects
    manually. You will usually get them from
    :meth:`~sage.libs.ppl.Polyhedron.relation_with`. You can also get
    pre-defined relations from the class methods :meth:`nothing`,
    :meth:`is_disjoint`, :meth:`strictly_intersects`,
    :meth:`is_included`, and :meth:`saturates`.

    Examples:

    >>> from ppl import Poly_Con_Relation
    >>> saturates     = Poly_Con_Relation.saturates();  saturates
    saturates
    >>> is_included   = Poly_Con_Relation.is_included(); is_included
    is_included
    >>> is_included.implies(saturates)
    False
    >>> saturates.implies(is_included)
    False
    >>> rels = []
    >>> rels.append(Poly_Con_Relation.nothing())
    >>> rels.append(Poly_Con_Relation.is_disjoint())
    >>> rels.append(Poly_Con_Relation.strictly_intersects())
    >>> rels.append(Poly_Con_Relation.is_included())
    >>> rels.append(Poly_Con_Relation.saturates())
    >>> rels
    [nothing, is_disjoint, strictly_intersects, is_included, saturates]
    >>> for i, rel_i in enumerate(rels):
    ...       for j, rel_j in enumerate(rels):
    ...           print int(rel_i.implies(rel_j)),
    ...       print
    1 0 0 0 0
    1 1 0 0 0
    1 0 1 0 0
    1 0 0 1 0
    1 0 0 0 1
    """
    def __cinit__(self, do_not_construct_manually=False):
        """
        The Cython constructor.

        See :class:`Poly_Con_Relation` for documentation.

        Tests:

        >>> from ppl import Poly_Con_Relation
        >>> Poly_Con_Relation.nothing()
        nothing
        """
        assert(do_not_construct_manually)
        self.thisptr = NULL

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        assert self.thisptr!=NULL, 'Do not construct Poly_Con_Relation objects manually!'
        del self.thisptr

    def implies(self, Poly_Con_Relation y):
        r"""
        Test whether ``self`` implies ``y``.

        INPUT:

        - ``y`` -- a :class:`Poly_Con_Relation`.

        OUTPUT:

        Boolean. ``True`` if and only if ``self`` implies ``y``.

        Examples:

            >>> from ppl import Poly_Con_Relation
            >>> nothing = Poly_Con_Relation.nothing()
            >>> nothing.implies( nothing )
            True
        """
        return self.thisptr.implies(y.thisptr[0])

    @classmethod
    def nothing(cls):
        r"""
        Return the assertion that says nothing.

        OUTPUT:

        A :class:`Poly_Con_Relation`.

        Examples:

        >>> from ppl import Poly_Con_Relation
        >>> Poly_Con_Relation.nothing()
        nothing
        """
        return _wrap_Poly_Con_Relation(PPL_Poly_Con_Relation_nothing())

    @classmethod
    def is_disjoint(cls):
        r"""
        Return the assertion "The polyhedron and the set of points
        satisfying the constraint are disjoint".

        OUTPUT:

        A :class:`Poly_Con_Relation`.

        Examples:

        >>> from ppl import Poly_Con_Relation
        >>> Poly_Con_Relation.is_disjoint()
        is_disjoint
        """
        return _wrap_Poly_Con_Relation(PPL_Poly_Con_Relation_is_disjoint())

    @classmethod
    def strictly_intersects(cls):
        r"""
        Return the assertion "The polyhedron intersects the set of
        points satisfying the constraint, but it is not included in
        it".

        :return: a :class:`Poly_Con_Relation`.

        Examples:

        >>> from ppl import Poly_Con_Relation
        >>> Poly_Con_Relation.strictly_intersects()
        strictly_intersects
        """
        return _wrap_Poly_Con_Relation(PPL_Poly_Con_Relation_strictly_intersects())


    @classmethod
    def is_included(cls):
        r"""
        Return the assertion "The polyhedron is included in the set of
        points satisfying the constraint".

        OUTPUT:

        A :class:`Poly_Con_Relation`.

        Examples:

        >>> from ppl import Poly_Con_Relation
        >>> Poly_Con_Relation.is_included()
        is_included
        """
        return _wrap_Poly_Con_Relation(PPL_Poly_Con_Relation_is_included())

    @classmethod
    def saturates(cls):
        r"""
        Return the assertion "".

        OUTPUT:

        A :class:`Poly_Con_Relation`.

        Examples:

        >>> from ppl import Poly_Con_Relation
        >>> Poly_Con_Relation.saturates()
        saturates
        """
        return _wrap_Poly_Con_Relation(PPL_Poly_Con_Relation_saturates())

    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        Examples:

        >>> cmd  = 'from ppl import Poly_Con_Relation\n'
        >>> cmd += 'Poly_Con_Relation.nothing().ascii_dump()\n'
        >>> import subprocess
        >>> proc = subprocess.Popen(['python', '-c', cmd], stderr=subprocess.PIPE)
        >>> out, err = proc.communicate()
        >>> print err
        NOTHING
        """
        self.thisptr.ascii_dump()

    def OK(self, check_non_empty=False):
        """
        Check if all the invariants are satisfied.

        Examples:

        >>> from ppl import Poly_Con_Relation
        >>> Poly_Con_Relation.nothing().OK()
        True
        """
        return self.thisptr.OK()

    def __repr__(self):
        r"""
        Return a string representation.

        OUTPUT:

        String.

        Examples:

        >>> from ppl import Poly_Con_Relation
        >>> Poly_Con_Relation.nothing().__repr__()
        'nothing'
        """
        rel = []
        if self.implies(Poly_Con_Relation.is_disjoint()):
            rel.append('is_disjoint')
        if self.implies(Poly_Con_Relation.strictly_intersects()):
            rel.append('strictly_intersects')
        if self.implies(Poly_Con_Relation.is_included()):
            rel.append('is_included')
        if self.implies(Poly_Con_Relation.saturates()):
            rel.append('saturates')
        if rel:
            return ', '.join(rel)
        else:
            return 'nothing'



