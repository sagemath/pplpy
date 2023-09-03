# distutils: language = c++
# distutils: libraries = gmp gmpxx ppl m
#*****************************************************************************
#       Copyright (C) 2010-2014 Volker Braun  <vbraun.name@gmail.com>
#                     2011 Simon King <simon.king@uni-jena.de>
#                     2011-2017 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#                     2012 Risan <ptrrsn.1@gmail.com>
#                     2013 Julien Puydt <julien.puydt@laposte.net>
#                     2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#                     2015 André Apitzsch <andre.apitzsch@etit.tu-chemnitz.de>
#                     2016 Jori Mäntysalo <jori.mantysalo@uta.fi>
#                     2016 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#                     2016-2017 Frédéric Chapoton <chapoton@math.univ-lyon1.fr>
#                     2016-2018 Vincent Delecroix <vincent.delecroix@labri.fr>
#                     2017-2018 Vincent Klein <vinklein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_NE, Py_GT, Py_GE
from cysignals.signals cimport sig_on, sig_off
from gmpy2 cimport GMPy_MPZ_From_mpz, import_gmpy2

from .ppl_decl cimport *
from .constraint cimport Constraint_System, Poly_Con_Relation, _wrap_Constraint_System
from .generator cimport Generator_System, Poly_Gen_Relation, _wrap_Generator_System
from .linear_algebra cimport *

# PPL can use floating-point arithmetic to compute integers
cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library":
    cdef void set_rounding_for_PPL()
    cdef void restore_pre_PPL_rounding()

# initialize gmpy2 C API
import_gmpy2()

# but with PPL's rounding the gsl will be very unhappy; must turn off!
restore_pre_PPL_rounding()


cdef class Polyhedron(object):
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

        >>> from ppl.polyhedron import Polyhedron
        >>> Polyhedron()
        Traceback (most recent call last):
        ...
        NotImplementedError: The Polyhedron class is abstract, you must not instantiate it.
        """
        raise NotImplementedError('The Polyhedron class is abstract, you must not instantiate it.')

    def __hash__(self):
        r"""
        Tests:

        >>> import ppl
        >>> hash(ppl.C_Polyhedron(ppl.point()))
        Traceback (most recent call last):
        ...
        TypeError: Polyhedron not hashable

        >>> cs = ppl.Constraint_System(ppl.Variable(0) > 0)
        >>> hash(NNC_Polyhedron(cs))
        Traceback (most recent call last):
        ...
        TypeError: Polyhedron not hashable
        """
        raise TypeError('Polyhedron not hashable')

    def __repr__(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        Examples:

        >>> from ppl import Variable, C_Polyhedron
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> C_Polyhedron( 5*x-2*y >=  x+y-1 ).__repr__()
        'A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 ray, 1 line'

        Special cases::

        >>> C_Polyhedron(3, 'empty').__repr__()
        'The empty polyhedron in QQ^3'
        >>> C_Polyhedron(3, 'universe').__repr__()
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
            if n_points == 1:
                desc += ' point'
            else:
                desc += ' points'
        if n_closure_points>0:
            if not first:
                desc += ", "
            first = False
            desc += str(n_closure_points)
            if n_closure_points == 1:
                desc += ' closure_point'
            else:
                desc += ' closure_points'
        if n_rays>0:
            if not first:
                desc += ", "
            first = False
            desc += str(n_rays)
            if n_rays == 1:
                desc += ' ray'
            else:
                desc += ' rays'
        if n_lines>0:
            if not first:
                desc += ", "
            first = False
            desc += repr(n_lines)
            if n_lines == 1:
                desc +=' line'
            else:
                desc +=' lines'
        return desc

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

        >>> from ppl import Variable, C_Polyhedron, point
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

        >>> from ppl import Poly_Con_Relation
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
        <BLANKLINE>
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

        >>> from ppl import Variable, C_Polyhedron, point
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

        >>> from ppl import Variable, C_Polyhedron, NNC_Polyhedron, Constraint_System
        >>> x = Variable(0);  y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert(x >= 0)
        >>> cs.insert(y >= 0)
        >>> cs.insert(3*x+5*y <= 10)
        >>> p = C_Polyhedron(cs)
        >>> pm = p.maximize(x+y)
        >>> for key in sorted(pm):
        ...     print("{} {}".format(key, pm[key]))
        bounded True
        generator point(10/3, 0/3)
        maximum True
        sup_d 3
        sup_n 10

        Unbounded case:

        >>> cs = Constraint_System()
        >>> cs.insert(x > 0)
        >>> p = NNC_Polyhedron(cs)
        >>> p.maximize(+x)
        {'bounded': False}
        >>> pm = p.maximize(-x)
        >>> for key in sorted(pm):
        ...     print("{} {}".format(key, pm[key]))
        bounded True
        generator closure_point(0/1)
        maximum False
        sup_d 1
        sup_n 0
        """
        cdef PPL_Coefficient sup_n
        cdef PPL_Coefficient sup_d
        cdef Generator g = Generator.point()
        cdef cppbool maximum
        sig_on()
        rc = self.thisptr.maximize(<PPL_Linear_Expression&>expr.thisptr[0], sup_n, sup_d, maximum, g.thisptr[0])
        sig_off()

        mpz_sup_n = GMPy_MPZ_From_mpz(sup_n.get_mpz_t())
        mpz_sup_d = GMPy_MPZ_From_mpz(sup_d.get_mpz_t())

        if rc:
            return {'bounded': True, 'sup_n': mpz_sup_n, 'sup_d': mpz_sup_d, 'maximum': maximum, 'generator': g}
        else:
            return {'bounded': False}

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

        >>> from ppl import Variable, C_Polyhedron, NNC_Polyhedron, Constraint_System
        >>> x = Variable(0);  y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x>=0 )
        >>> cs.insert( y>=0 )
        >>> cs.insert( 3*x+5*y<=10 )
        >>> p = C_Polyhedron(cs)
        >>> pm = p.minimize( x+y )
        >>> for key in sorted(pm):
        ...     print("{} {}".format(key, pm[key]))
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
        ...    print("{} {}".format(key, pm[key]))
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
        cdef Generator g = Generator.point()
        cdef cppbool minimum
        sig_on()
        rc = self.thisptr.minimize(<PPL_Linear_Expression&>expr.thisptr[0], inf_n, inf_d, minimum, g.thisptr[0])
        sig_off()

        mpz_inf_n = GMPy_MPZ_From_mpz(inf_n.get_mpz_t())
        mpz_inf_d = GMPy_MPZ_From_mpz(inf_d.get_mpz_t())

        if rc:
            return {'bounded': True, 'inf_n': mpz_inf_n, 'inf_d': mpz_inf_d, 'minimum': minimum, 'generator': g}
        else:
            return {'bounded': False}

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
        <BLANKLINE>
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
        <BLANKLINE>
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

        See also :meth:`add_constraints`.

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

        See also :meth:`add_constraint`.

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

        See also :meth:`add_generator`.

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
        <BLANKLINE>
        """
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
        <BLANKLINE>
        """
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

        >>> from ppl import Variable, C_Polyhedron, point, NNC_Polyhedron
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
        <BLANKLINE>
        """
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

            The modified polyhedron is not necessarily a lattice
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
        sig_on()
        self.thisptr.topological_closure_assign()
        sig_off()

    def BHRZ03_widening_assign(self, Polyhedron y, unsigned tp = 0):
        r"""
        Assigns to ``self``` the result of computing the `BHRZ03-widening`
        between ``self`` and ``y``.

        INPUT:

        - ``y`` -- a :class:`Polyhedron` that must be contained in ``self``

        - ``tp`` -- an optional unsigned variable with the number of
          available tokens (to be used when applying the `widening with
          tokens` delay technique).

        OUTPUT:

        This method assigns to ``self`` the result of computing the
        BHRZ03-widening between ``self`` and ``y``. And returns the
        new value of ``tp``.

        Raises a ``ValueError`` if ``self`` and ``y`` are
        topology-incompatible.

        Examples:

        >>> from ppl import NNC_Polyhedron, Variable
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> ph1 = NNC_Polyhedron(2)
        >>> ph1.add_constraint( y >= 0 )
        >>> ph1.add_constraint( x + y > 0 )
        >>> ph1.add_constraint( x - y < 1 )
        >>> ph2 = NNC_Polyhedron(2)
        >>> ph2.add_constraint( y >= 0 )
        >>> ph2.add_constraint( x > 0 )
        >>> ph2.add_constraint( x < 1 )
        >>> tp = ph1.BHRZ03_widening_assign( ph2 )
        >>> known_result = NNC_Polyhedron(2)
        >>> known_result.add_constraint(y >= 0)
        >>> known_result == ph1
        True

        >>> from ppl import C_Polyhedron, Generator_System, point, ray
        >>> gs1 = Generator_System()
        >>> gs1.insert(point())
        >>> gs1.insert(point( x + 2 * y ))
        >>> gs1.insert(ray( x ))
        >>> gs1.insert(ray( 2 * x + y ))
        >>> ph1 = C_Polyhedron(gs1)
        >>> gs2 = Generator_System()
        >>> gs2.insert(point())
        >>> gs2.insert(point( x + 2 * y ))
        >>> gs2.insert(ray( x ))
        >>> gs2.insert(ray( x + y ))
        >>> ph2 = C_Polyhedron( gs2 )
        >>> tp = ph2.BHRZ03_widening_assign( ph1 )
        >>> known_result = C_Polyhedron(2)
        >>> known_result.add_constraint( y >= 0 )
        >>> known_result.add_constraint( 2 * x - y >= 0 )
        >>> ph2 == known_result
        True

        if this method is going to be called too many times, the
        ``tp`` parameter allows one to reduce the computation cost by
        computing the widening only when ``tp`` is equal to 0,
        otherwise the this method will decrement the value of ``tp``
        and return it.

        >>> from ppl import closure_point
        >>> gs1 = Generator_System()
        >>> gs1.insert(point( 2 * x ))
        >>> gs1.insert(closure_point( x + y ))
        >>> gs1.insert(closure_point( 3 * x + y ))
        >>> ph1 = NNC_Polyhedron(gs1)
        >>> gs2 = Generator_System()
        >>> gs2.insert(point( 2 * x ))
        >>> gs2.insert(closure_point( y ))
        >>> gs2.insert(closure_point( 4 * x + y ))
        >>> ph2 = NNC_Polyhedron(gs2)
        >>> ph2_copy = NNC_Polyhedron(ph2)
        >>> known_result = NNC_Polyhedron(2)
        >>> known_result.add_constraint(y >= 0)
        >>> known_result.add_constraint(y < 1)
        >>> tp = ph2.BHRZ03_widening_assign(ph1, 1)
        >>> tp == 0
        True
        >>> ph2 == ph2_copy
        True
        >>> tp = ph2.BHRZ03_widening_assign(ph1, 0)
        >>> tp == 0
        True
        >>> ph2 == known_result
        True
        """
        sig_on()
        try:
            self.thisptr.BHRZ03_widening_assign(y.thisptr[0], &tp)
        finally:
            sig_off()
        return tp

    def limited_BHRZ03_extrapolation_assign(self, Polyhedron y, Constraint_System cs, unsigned tp = 0):
        r"""
        Assigns to ``self``` the result of computing the limited
        extrapolation between ``self`` and ``y`` using the
        `BHRZ03-widening` operator.

        INPUT:

        - ``y`` -- a :class:`Polyhedron` that must be contained in ``self``

        - ``cs`` -- a :class:`Constraint_System` used to improve the
          widened polyhedron

        - ``tp`` -- an optional unsigned variable with the number of
          available tokens (to be used when applying the `widening with
          tokens` delay technique).

        OUTPUT:

        This method assigns to ``self`` the result of computing the
        limited extrapolation between ``self`` and ``y`` using the
        BHRZ03-widening operator. And returns the new value of ``tp``.

        Raises a ``ValueError`` if ``self`` and ``y`` are
        topology-incompatible or dimesion-incompatible.

        Examples:
        >>> from ppl import C_Polyhedron, Generator_System, Variable, point, Constraint_System
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> gs1 = Generator_System()
        >>> gs1.insert(point())
        >>> gs1.insert(point( x + y ))
        >>> gs1.insert(point( x ))
        >>> ph1 = C_Polyhedron( gs1 )
        >>> gs2 = Generator_System()
        >>> gs2.insert(point())
        >>> gs2.insert(point( 2 * x ))
        >>> gs2.insert(point( 2 * x + 2 * y ))
        >>> ph2 = C_Polyhedron( gs2 )
        >>> cs = Constraint_System()
        >>> cs.insert( x <= 5 )
        >>> cs.insert( y <= 4 )
        >>> tp = ph2.limited_BHRZ03_extrapolation_assign(ph1, cs)
        >>> known_result = C_Polyhedron(2)
        >>> known_result.add_constraint(y >= 0)
        >>> known_result.add_constraint(x - y >= 0)
        >>> known_result.add_constraint(y <= 4)
        >>> known_result.add_constraint(x <= 5)
        >>> known_result == ph2
        True
        """
        sig_on()
        try:
            self.thisptr.limited_BHRZ03_extrapolation_assign(y.thisptr[0],
                                                             cs.thisptr[0],
                                                             &tp)
        finally:
            sig_off()
        return tp

    def bounded_BHRZ03_extrapolation_assign(self, Polyhedron y, Constraint_System cs, unsigned tp = 0):
        r"""
        Assigns to ``self``` the result of computing the bounded
        extrapolation between ``self`` and ``y`` using the
        `BHRZ03-widening` operator.

        INPUT:

        - ``y`` -- a :class:`Polyhedron` that must be contained in ``self``

        - ``cs`` -- a :class:`Constraint_System` used to improve the
          widened polyhedron

        - ``tp`` -- an optional unsigned variable with the number of
          available tokens (to be used when applying the `widening with
          tokens` delay technique).

        OUTPUT:

        This method assigns to ``self`` the result of computing the
        bounded extrapolation between ``self`` and ``y`` using the
        BHRZ03-widening operator. And returns the new value of ``tp``.

        Raises a ``ValueError`` if ``self`` and ``y`` are
        topology-incompatible or dimesion-incompatible.

        Examples:

        >>> from ppl import Variable, Constraint_System, C_Polyhedron
        >>> x = Variable(0)
        >>> ph1 = C_Polyhedron(1)
        >>> ph1.add_constraint( 1 <= x )
        >>> ph1.add_constraint( x <= 2 )
        >>> ph2 = C_Polyhedron(1)
        >>> ph2.add_constraint( 0 <= x )
        >>> ph2.add_constraint( x <= 3 )
        >>> cs = Constraint_System()
        >>> tp = ph2.bounded_BHRZ03_extrapolation_assign(ph1, cs)
        >>> known_result = C_Polyhedron(1)
        >>> known_result.add_constraint(0 <= x)
        >>> known_result == ph2
        True
        """
        sig_on()
        try:
            self.thisptr.bounded_BHRZ03_extrapolation_assign(y.thisptr[0],
                                                             cs.thisptr[0],
                                                             &tp)
        finally:
            sig_off()
        return tp

    def H79_widening_assign(self, Polyhedron y, unsigned tp = 0):
        r"""
        Assigns to ``self``` the result of computing the `H79-widening`
        between ``self`` and ``y``.

        INPUT:

        - ``y`` -- a :class:`Polyhedron` that must be contained in ``self``

        - ``tp`` -- an optional unsigned variable with the number of
          available tokens (to be used when applying the `widening with
          tokens` delay technique).

        OUTPUT:

        This method assigns to ``self`` the result of computing the
        BHRZ03-widening between ``self`` and ``y``. And returns the
        new value of ``tp``.

        Raises a ``ValueError`` if ``self`` and ``y`` are
        topology-incompatible.


        Examples:
        >>> from ppl import Variable, C_Polyhedron
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> ph1 = C_Polyhedron(2)
        >>> ph1.add_constraint( x >= 2 )
        >>> ph1.add_constraint( y >= 0 )
        >>> ph2 = C_Polyhedron(2)
        >>> ph2.add_constraint( x >= 0 )
        >>> ph2.add_constraint( y >= 0 )
        >>> ph2.add_constraint( x-y >= 2 )
        >>> tp = ph1.H79_widening_assign(ph2)
        >>> known_result = C_Polyhedron(2)
        >>> known_result.add_constraint(y >= 0)
        >>> known_result == ph1
        True

        ``self`` and ``y`` must be dimension- and topology-compatible,
        or an exception is raised:

        >>> z = Variable(2)
        >>> ph1.H79_widening_assign( C_Polyhedron(z>=0) )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::H79_widening_assign(y):
        this->space_dimension() == 2, y.space_dimension() == 3.

        if this method is going to be called too many times, the
        ``tp`` parameter allows one to reduce the computation cost by
        computing the widening only when ``tp`` is equal to 0,
        otherwise the this method will decrement the value of ``tp``
        and return it.

        >>> from ppl import point, ray, Generator_System
        >>> gs1 = Generator_System()
        >>> gs1.insert(point())
        >>> gs1.insert(ray( x + y ))
        >>> gs1.insert(ray( x ))
        >>> ph1 = C_Polyhedron(gs1)
        >>> gs2 = Generator_System()
        >>> gs2.insert(point())
        >>> gs2.insert(ray( x ))
        >>> gs2.insert(ray( x + 2 * y ))
        >>> ph2 = C_Polyhedron(gs2)
        >>> ph2_copy = C_Polyhedron(ph2)
        >>> known_result = C_Polyhedron(2)
        >>> known_result.add_constraint(y >= 0)
        >>> tp = ph2.H79_widening_assign(ph1, 1)
        >>> tp == 0
        True
        >>> ph2 == ph2_copy
        True
        >>> tp = ph2.H79_widening_assign(ph1, 0)
        >>> tp == 0
        True
        >>> ph2 == known_result
        True
        """
        sig_on()
        try:
            self.thisptr.H79_widening_assign(y.thisptr[0], &tp)
        finally:
            sig_off()
        return tp

    widening_assign = H79_widening_assign

    def limited_H79_extrapolation_assign(self, Polyhedron y, Constraint_System cs, unsigned tp = 0):
        r"""
        Assigns to ``self``` the result of computing the limited
        extrapolation between ``self`` and ``y`` using the
        `H79-widening` operator.

        INPUT:

        - ``y`` -- a :class:`Polyhedron` that must be contained in ``self``

        - ``cs`` -- a :class:`Constraint_System` used to improve the
          widened polyhedron

        - ``tp`` -- an optional unsigned variable with the number of
          available tokens (to be used when applying the `widening with
          tokens` delay technique).

        OUTPUT:

        This method assigns to ``self`` the result of computing the
        limited extrapolation between ``self`` and ``y`` using the
        H79-widening operator. And returns the new value of ``tp``.

        Raises a ``ValueError`` if ``self`` and ``y`` are
        topology-incompatible or dimesion-incompatible.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, Constraint_System, point
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> ph1 = C_Polyhedron(2)
        >>> ph1.add_constraint( x >= 0 )
        >>> ph1.add_constraint( x <= 1 )
        >>> ph1.add_constraint( y >= 0 )
        >>> ph1.add_constraint( x - y >= 0 )
        >>> ph2 = C_Polyhedron(2)
        >>> ph2.add_constraint( x >= 0 )
        >>> ph2.add_constraint( x <= 2 )
        >>> ph2.add_constraint( y >= 0 )
        >>> ph2.add_constraint( x - y >= 0 )
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0 )
        >>> cs.insert( y >= 0 )
        >>> cs.insert( x <= 5 )
        >>> cs.insert( y <= 5 )
        >>> tp = ph2.limited_H79_extrapolation_assign(ph1, cs)
        >>> known_result = C_Polyhedron(2)
        >>> known_result.add_constraint( x - y >= 0)
        >>> known_result.add_constraint(y >= 0)
        >>> known_result.add_constraint(x <= 5)
        >>> known_result == ph2
        True

        >>> ph1 = C_Polyhedron(2)
        >>> ph1.add_constraint( x >= 0 )
        >>> ph1.add_constraint( x <= 1 )
        >>> ph1.add_constraint( y == 0 )
        >>> ph2 = C_Polyhedron(2)
        >>> ph2.add_constraint( x <= 2 )
        >>> ph2.add_constraint( y >= 0 )
        >>> ph2.add_constraint( y <= x )
        >>> cs = Constraint_System()
        >>> cs.insert( x <= 5 )
        >>> cs.insert( y <= -1 )
        >>> known_result = C_Polyhedron(ph2)
        >>> known_result.add_generator(point(5*x))
        >>> known_result.add_generator(point(5*x + 5*y))
        >>> tp = ph2.limited_H79_extrapolation_assign(ph1, cs)
        >>> known_result == ph2
        True
        """
        sig_on()
        try:
            self.thisptr.limited_H79_extrapolation_assign(y.thisptr[0],
                                                          cs.thisptr[0],
                                                          &tp)
        finally:
            sig_off()
        return tp

    def bounded_H79_extrapolation_assign(self, Polyhedron y, Constraint_System cs, unsigned tp = 0):
        r"""
        Assigns to ``self``` the result of computing the bounded
        extrapolation between ``self`` and ``y`` using the
        `H79-widening` operator.

        INPUT:

        - ``y`` -- a :class:`Polyhedron` that must be contained in ``self``

        - ``cs`` -- a :class:`Constraint_System` used to improve the
          widened polyhedron

        - ``tp`` -- an optional unsigned variable with the number of
          available tokens (to be used when applying the `widening with
          tokens` delay technique).

        OUTPUT:

        This method assigns to ``self`` the result of computing the
        bounded extrapolation between ``self`` and ``y`` using the
        H79-widening operator. And returns the new value of ``tp``.

        Raises a ``ValueError`` if ``self`` and ``y`` are
        topology-incompatible or dimesion-incompatible.

        Examples:

        >>> from ppl import Variable, C_Polyhedron, Constraint_System
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> ph1 = C_Polyhedron(2)
        >>> ph1.add_constraint( x-3 >= 0 )
        >>> ph1.add_constraint( x-3 <= 1 )
        >>> ph1.add_constraint( y >= 0 )
        >>> ph1.add_constraint( y <= 1 )
        >>> ph2 = C_Polyhedron(2)
        >>> ph2.add_constraint( 2*x-5 >= 0 )
        >>> ph2.add_constraint( x-3 <= 1 )
        >>> ph2.add_constraint( 2*y+3 >= 0 )
        >>> ph2.add_constraint( 2*y-5 <= 0 )
        >>> cs = Constraint_System()
        >>> cs.insert( x >= y )
        >>> tp = ph2.bounded_H79_extrapolation_assign(ph1, cs)
        >>> known_result = C_Polyhedron(2)
        >>> known_result.add_constraint( x >= 2 )
        >>> known_result.add_constraint( x <= 4 )
        >>> known_result.add_constraint( y >= -2 )
        >>> known_result.add_constraint( x >= y )
        >>> known_result == ph2
        True
        """
        sig_on()
        try:
            self.thisptr.bounded_H79_extrapolation_assign(y.thisptr[0],
                                                          cs.thisptr[0],
                                                          &tp)
        finally:
            sig_off()
        return tp

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
            (x,y,z)^T \in \mathbb{R}^3
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
            (x,y,0)^T \in \mathbb{R}^3
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
        m = int(m)
        sig_on()
        try:
            self.thisptr.add_space_dimensions_and_project(m)
        finally:
            sig_off()

    def concatenate_assign(self, Polyhedron y):
        r"""
        Assign to ``self`` the concatenation of ``self`` and ``y``.

        This functions returns the Cartesian product of ``self`` and
        ``y``.

        Viewing a polyhedron as a set of tuples (its points), it is
        sometimes useful to consider the set of tuples obtained by
        concatenating an ordered pair of polyhedra. Formally, the
        concatenation of the polyhedra `P` and `Q` (taken in this
        order) is the polyhedron such that

        .. MATH::

            R =
            \Big\{
            (x_0,\dots,x_{n-1},y_0,\dots,y_{m-1})^T \in \mathbb{R}^{n+m}
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
        <BLANKLINE>
        >>> p1.concatenate_assign( C_Polyhedron(p1.max_space_dimension(), 'empty') )
        Traceback (most recent call last):
        ...
        ValueError: PPL::C_Polyhedron::concatenate_assign(y):
        concatenation exceeds the maximum allowed space dimension.
        """
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
        new_dimension = int(new_dimension)
        sig_on()
        try:
            self.thisptr.remove_higher_space_dimensions(new_dimension)
        finally:
            sig_off()

    def affine_image(self, Variable v, Linear_Expression le):
        r"""
        Set this polyhedron to the image of the map `v -> le`

        INPUT:

        - ``v`` -- a variable

        - ``le`` -- a linear expression

        Examples:

        >>> import ppl
        >>> x = ppl.Variable(0)
        >>> y = ppl.Variable(1)
        >>> gs0 = ppl.Generator_System()
        >>> gs0.insert(ppl.point())
        >>> gs0.insert(ppl.point(x))
        >>> gs0.insert(ppl.point(y))
        >>> gs0.insert(ppl.point(x+y))
        >>> p0 = ppl.C_Polyhedron(gs0)
        >>> gs1 = ppl.Generator_System()
        >>> gs1.insert(ppl.point())
        >>> gs1.insert(ppl.point(x))
        >>> gs1.insert(ppl.point(x+y))
        >>> gs1.insert(ppl.point(2*x+y))
        >>> p1 = ppl.C_Polyhedron(gs1)

        >>> p0 == p1
        False
        >>> p0.affine_image(x, x+y)
        >>> p0 == p1
        True
        """
        self.thisptr.affine_image(v.thisptr[0], le.thisptr[0])

    def affine_preimage(self, Variable v, Linear_Expression le):
        r"""
        Set this polyhedron to the preimage of the map `v -> le`

        INPUT:

        - ``v`` -- a variable

        - ``le`` -- a linear expression

        Examples:

        >>> import ppl
        >>> x = ppl.Variable(0)
        >>> y = ppl.Variable(1)
        >>> gs0 = ppl.Generator_System()
        >>> gs0.insert(ppl.point())
        >>> gs0.insert(ppl.point(x))
        >>> gs0.insert(ppl.point(y))
        >>> gs0.insert(ppl.point(x+y))
        >>> p0 = ppl.C_Polyhedron(gs0)
        >>> gs1 = ppl.Generator_System()
        >>> gs1.insert(ppl.point())
        >>> gs1.insert(ppl.point(x))
        >>> gs1.insert(ppl.point(x+y))
        >>> gs1.insert(ppl.point(2*x+y))
        >>> p1 = ppl.C_Polyhedron(gs1)

        >>> p0 == p1
        False
        >>> p1.affine_preimage(x, x+y)
        >>> p0 == p1
        True
        """
        self.thisptr.affine_preimage(v.thisptr[0], le.thisptr[0])

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
        >>> from subprocess import Popen, PIPE
        >>> import sys
        >>> proc = Popen([sys.executable, '-c', cmd], stdout=PIPE, stderr=PIPE)
        >>> out, err = proc.communicate()
        >>> len(out)
        0
        >>> print(str(err.decode('ascii')))
        space_dim 2
        ...
        con_sys (up-to-date)
        topology NECESSARILY_CLOSED
        ...
        <BLANKLINE>
        sat_c
        0 x 0
        <BLANKLINE>
        sat_g
        2 x 2
        0 0
        0 1
        <BLANKLINE>
        <BLANKLINE>
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
        >>> C_Polyhedron(1, 'empty').max_space_dimension() > 2**20
        True
        """
        return self.thisptr.max_space_dimension()

    def hash_code(self):
        r"""
        Return a hash code

        Tests:

        >>> from ppl import Constraint_System, Variable, C_Polyhedron
        >>> x = Variable(0)
        >>> p = C_Polyhedron( 5*x >= 3 )
        >>> p.hash_code()
        1

        >>> y = Variable(1)
        >>> cs = Constraint_System()
        >>> cs.insert( x >= 0 )
        >>> cs.insert( y >= 0 )
        >>> p = C_Polyhedron(cs)
        >>> p.hash_code()
        2
        """
        return self.thisptr[0].hash_code()

    def __richcmp__(Polyhedron lhs, Polyhedron rhs, int op):
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
        if op == Py_LT:      # <   0
            result = rhs.strictly_contains(lhs)
        elif op == Py_LE:    # <=  1
            result = rhs.contains(lhs)
        elif op == Py_EQ:    # ==  2
            result = (lhs.thisptr[0] == rhs.thisptr[0])
        elif op == Py_GT:    # >   4
            result = lhs.strictly_contains(rhs)
        elif op == Py_GE:    # >=  5
            result = lhs.contains(rhs)
        elif op == Py_NE:    # !=  3
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

    >>> from ppl import Constraint_System, Generator_System, Variable, C_Polyhedron, point, ray
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
            raise ValueError('Cannot initialize C_Polyhedron with '+str(arg)+'.')
        degenerate_element = degenerate_element.lower()
        if degenerate_element=='universe':
            self.thisptr = new PPL_C_Polyhedron(<PPL_dimension_type>dim, UNIVERSE)
            return
        elif degenerate_element=='empty':
            self.thisptr = new PPL_C_Polyhedron(<PPL_dimension_type>dim, EMPTY)
            return
        else:
            raise ValueError('Unknown value: degenerate_element='+str(degenerate_element)+'.')

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
        >>> from pickle import loads, dumps
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


####################################################
### NNC_Polyhedron ###################################
####################################################
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
        >>> from pickle import loads, dumps
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
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 closure_point, 1 ray, 1 line
        """
        if self.is_empty():
            return (NNC_Polyhedron, (self.space_dimension(), 'empty'))
        elif self.is_universe():
            return (NNC_Polyhedron, (self.space_dimension(), 'universe'))
        else:
            return (NNC_Polyhedron, (self.generators(),))
