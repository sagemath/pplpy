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

from gmpy2 cimport GMPy_MPZ_From_mpz, import_gmpy2
from cython.operator cimport dereference as deref
from .linear_algebra cimport PPL_Coefficient_from_pyobject


# PPL can use floating-point arithmetic to compute integers
cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library":
    cdef void set_rounding_for_PPL()
    cdef void restore_pre_PPL_rounding()

# initialize gmpy2 C API
import_gmpy2()

# but with PPL's rounding the gsl will be very unhappy; must turn off!
restore_pre_PPL_rounding()

####################################################
cdef class Generator(object):
    r"""
    Wrapper for PPL's ``Generator`` class.

    An object of the class Generator is one of the following:

      * a line :math:`\ell = (a_0, \dots, a_{n-1})^T`

      * a ray :math:`r = (a_0, \dots, a_{n-1})^T`

      * a point :math:`p = (\tfrac{a_0}{d}, \dots, \tfrac{a_{n-1}}{d})^T`

      * a closure point :math:`c = (\tfrac{a_0}{d}, \dots, \tfrac{a_{n-1}}{d})^T`

    where :math:`n` is the dimension of the space and, for points and
    closure points, :math:`d` is the divisor.

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

    def __hash__(self):
        r"""
        Tests:

        >>> import ppl
        >>> hash(ppl.point())
        Traceback (most recent call last):
        ...
        TypeError: Generator unhashable
        """
        raise TypeError('Generator unhashable')

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
            g.thisptr = new PPL_Generator(PPL_line(e.thisptr[0]))
        except BaseException:
            # g.thisptr must be set to something valid or g.__dealloc__() will segfault
            g.thisptr = new PPL_Generator(PPL_point(e.thisptr[0], PPL_Coefficient(1)))
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
            g.thisptr = new PPL_Generator(PPL_ray(e.thisptr[0]))
        except BaseException:
            # g.thisptr must be set to something valid or g.__dealloc__() will segfault
            g.thisptr = new PPL_Generator(PPL_point(e.thisptr[0], PPL_Coefficient(1)))
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
        # This does not work as Cython gets confused by the private default ctor
        #   return _wrap_Generator(PPL_point(e.thisptr[0], PPL_Coefficient(d.value)))
        # workaround follows
        cdef Generator g = Generator(True)
        try:
            g.thisptr = new PPL_Generator(PPL_point(e.thisptr[0], PPL_Coefficient_from_pyobject(divisor)))
        except BaseException:
            # g.thisptr must be set to something valid or g.__dealloc__() will segfault
            g.thisptr = new PPL_Generator(PPL_point(e.thisptr[0], PPL_Coefficient(1)))
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
        # This does not work as Cython gets confused by the private default ctor
        #   return _wrap_Generator(PPL_closure_point(e.thisptr[0], PPL_Coefficient(d.value)))
        # workaround follows
        cdef Generator g = Generator(True)
        try:
            g.thisptr = new PPL_Generator(PPL_closure_point(e.thisptr[0], PPL_Coefficient_from_pyobject(divisor)))
        except BaseException:
            # g.thisptr must be set to something valid or g.__dealloc__() will segfault
            g.thisptr = new PPL_Generator(PPL_point(e.thisptr[0], PPL_Coefficient(1)))
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

    def set_space_dimension(self, size_t n):
        r"""
        Set the dimension of this generator to ``n``

        Examples:

        >>> import ppl
        >>> p = ppl.point()
        >>> p
        point()
        >>> p.set_space_dimension(5)
        >>> p
        point(0/1, 0/1, 0/1, 0/1, 0/1)

        >>> p = ppl.point(1 * ppl.Variable(0) + 2 * ppl.Variable(1) + 3 * ppl.Variable(2))
        >>> p.set_space_dimension(2)
        >>> p
        point(1/1, 2/1)
        """
        self.thisptr.set_space_dimension(n)

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
        mpz(1)
        """
        return GMPy_MPZ_From_mpz(self.thisptr.coefficient(v.thisptr[0]).get_mpz_t())

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
            (mpz(3), mpz(5))
        """
        cdef int d = self.space_dimension()
        cdef int i
        coeffs = []
        for i in range(0, d):
            coeffs.append(GMPy_MPZ_From_mpz(self.thisptr.coefficient(PPL_Variable(i)).get_mpz_t()))
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
        mpz(1)
        >>> line = Generator.line(2*x-y+5)
        >>> line.divisor()
        Traceback (most recent call last):
        ...
        ValueError: PPL::Generator::divisor():
        *this is neither a point nor a closure point.
        """
        return GMPy_MPZ_From_mpz(self.thisptr.divisor().get_mpz_t())

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

        >>> from ppl import Variable, point, line
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
        >>> import sys
        >>> proc = subprocess.Popen([sys.executable, '-c', cmd], stderr=subprocess.PIPE)
        >>> out, err = proc.communicate()
        >>> print(str(err.decode('ascii')))
        size 3 1 3 2 ...
        <BLANKLINE>
        """
        self.thisptr.ascii_dump()

    def __reduce__(self):
        """
        Pickle object.

        Tests:

        >>> from ppl import Generator, Variable
        >>> from pickle import loads, dumps
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
        le = Linear_Expression(self.coefficients(), 0)
        t = self.thisptr.type()
        if t == LINE:
            return (Generator.line, (le,))
        elif t == RAY:
            return (Generator.ray, (le,))
        elif t == POINT:
            return (Generator.point, (le, self.divisor()))
        elif t == CLOSURE_POINT:
            return (Generator.closure_point, (le, self.divisor()))
        else:
            raise RuntimeError

    def permute_space_dimensions(self, cycle):
        r"""
        Permute the coordinates according to ``cycle``.

        >>> from ppl import Generator, Variable
        >>> x = Variable(0); y = Variable(1); z = Variable(2)
        >>> p = Generator.point(2*x+7*y-z, 3)
        >>> p.permute_space_dimensions([0, 1])
        >>> p
        point(7/3, 2/3, -1/3)
        >>> p.permute_space_dimensions([0, 2, 1])
        >>> p
        point(2/3, -1/3, 7/3)
        """
        cdef cppvector[PPL_Variable] cpp_cycle
        cdef int i
        for i in cycle:
            cpp_cycle.push_back(PPL_Variable(i))
        self.thisptr.permute_space_dimensions(cpp_cycle)

####################################################
### Generator_System  ##############################
####################################################
####################################################
cdef class Generator_System(object):
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
        if isinstance(arg, (list, tuple)):
            self.thisptr = new PPL_Generator_System()
            for generator in arg:
                self.insert(generator)
            return
        raise ValueError('Cannot initialize with '+str(arg)+'.')

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr

    def __hash__(self):
        r"""
        Tests:

        >>> import ppl
        >>> hash(ppl.Generator_System())
        Traceback (most recent call last):
        ...
        TypeError: Generator_System unhashable
        """
        raise TypeError('Generator_System unhashable')

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

    def set_space_dimension(self, size_t n):
        r"""
        Set the dimension of the vector space enclosing ``self``

        Examples:

        >>> import ppl
        >>> gs = ppl.Generator_System()
        >>> gs.insert(ppl.point(ppl.Variable(0) + ppl.Variable(12)))
        >>> gs.insert(ppl.point(ppl.Variable(1) + 2*ppl.Variable(3)))
        >>> gs.space_dimension()
        13
        >>> gs.set_space_dimension(25)
        >>> gs.space_dimension()
        25
        >>> gs.set_space_dimension(3)
        >>> gs.space_dimension()
        3
        """
        self.thisptr.set_space_dimension(n)

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
        >>> import subprocess, sys
        >>> proc = subprocess.Popen([sys.executable, '-c', cmd], stderr=subprocess.PIPE)
        >>> out, err = proc.communicate()
        >>> print(str(err.decode('ascii')))
        topology NECESSARILY_CLOSED
        ...
        <BLANKLINE>
        """
        self.thisptr.ascii_dump()

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

        >>> from ppl import Generator_System, Variable, point, ray, line, closure_point
        >>> x = Variable(0)
        >>> gs = Generator_System(point(3*x))
        >>> iter = gs.__iter__()
        >>> next(iter)
        point(3/1)
        >>> next(iter)
        Traceback (most recent call last):
        ...
        StopIteration
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> gs = Generator_System( line(5*x-2*y) )
        >>> gs.insert( ray(6*x-3*y) )
        >>> gs.insert( point(2*x-7*y, 5) )
        >>> gs.insert( closure_point(9*x-1*y, 2) )
        >>> next(gs.__iter__())
        line(5, -2)
        >>> list(gs)
        [line(5, -2), ray(2, -1), point(2/5, -7/5), closure_point(9/2, -1/2)]
        """
        cdef PPL_gs_iterator *gsi_ptr = new PPL_gs_iterator(self.thisptr[0].begin())
        try:
            while gsi_ptr[0] != self.thisptr[0].end():
                yield _wrap_Generator(deref(gsi_ptr[0].inc(1)))
        finally:
            del gsi_ptr

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
        s += ', '.join([repr(g) for g in self])
        s += '}'
        return s

    def __reduce__(self):
        """
        Pickle object.

        Tests:

        >>> from ppl import Generator_System, Variable, point, ray
        >>> from pickle import loads, dumps
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> gs = Generator_System((point(3*x+2*y+1), ray(x)));  gs
        Generator_System {point(3/1, 2/1), ray(1, 0)}
        >>> loads(dumps(gs))
        Generator_System {point(3/1, 2/1), ray(1, 0)}
        """
        return (Generator_System, (tuple(self), ))

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

cdef _wrap_Generator(PPL_Generator generator):
    """
    Wrap a C++ ``PPL_Generator`` into a Cython ``Generator``.
    """
    cdef Generator g = Generator(True)
    g.thisptr = new PPL_Generator(generator)
    return g


cdef _wrap_Generator_System(PPL_Generator_System generator_system):
    """
    Wrap a C++ ``PPL_Generator_System`` into a Cython ``Generator_System``.
    """
    cdef Generator_System gs = Generator_System()
    del gs.thisptr
    gs.thisptr = new PPL_Generator_System(generator_system)
    return gs

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

    def __hash__(self):
        r"""
        Tests:

        >>> import ppl
        >>> hash(ppl.Poly_Gen_Relation.nothing())
        Traceback (most recent call last):
        ...
        TypeError: Poly_Gen_Relation unhashable
        """
        raise TypeError('Poly_Gen_Relation unhashable')

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
        >>> from subprocess import Popen, PIPE
        >>> import sys
        >>> proc = Popen([sys.executable, '-c', cmd], stderr=PIPE)
        >>> out, err = proc.communicate()
        >>> print(str(err.decode('ascii')))
        NOTHING
        >>> proc.returncode
        0
        """
        self.thisptr.ascii_dump()

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

cdef _wrap_Poly_Gen_Relation(PPL_Poly_Gen_Relation relation):
    """
    Wrap a C++ ``PPL_Poly_Gen_Relation`` into a Cython ``Poly_Gen_Relation``.
    """
    cdef Poly_Gen_Relation rel = Poly_Gen_Relation(True)
    rel.thisptr = new PPL_Poly_Gen_Relation(relation)
    return rel
