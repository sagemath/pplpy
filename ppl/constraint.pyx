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
#                     2016-2020 Vincent Delecroix <vincent.delecroix@labri.fr>
#                     2017-2018 Vincent Klein <vinklein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from gmpy2 cimport GMPy_MPZ_From_mpz, import_gmpy2
from cython.operator cimport dereference as deref
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_NE, Py_GT, Py_GE

from .congruence cimport Congruence

# PPL can use floating-point arithmetic to compute integers
cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library":
    cdef void set_rounding_for_PPL()
    cdef void restore_pre_PPL_rounding()

# initialize gmpy2 C API
import_gmpy2()

# but with PPL's rounding the gsl will be very unhappy; must turn off!
restore_pre_PPL_rounding()


def _dummy():
    raise ValueError


cdef class Constraint(object):
    r"""
    Wrapper for PPL's ``Constraint`` class.

    An object of the class ``Constraint`` is either:

    * an equality :math:`\sum_{i=0}^{n-1} a_i x_i + b = 0`

    * a non-strict inequality :math:`\sum_{i=0}^{n-1} a_i x_i + b \geq 0`

    * a strict inequality :math:`\sum_{i=0}^{n-1} a_i x_i + b > 0`

    where :math:`n` is the dimension of the space, :math:`a_i` is the integer
    coefficient of variable :math:`x_i`, and :math:`b_i` is the integer
    inhomogeneous term.

    INPUT/OUTPUT:

    You construct constraints by writing inequalities in
    :class:`Linear_Expression`. Do not attempt to manually construct
    constraints.

    Examples:

    >>> from ppl import Variable, Linear_Expression
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
    mpz(-1)==0
    """
    def __init__(self, arg=None):
        if arg is None:
            self.thisptr = new PPL_Constraint()
        elif isinstance(arg, Constraint):
            self.thisptr = new PPL_Constraint((<Constraint> arg).thisptr[0])
        else:
            raise TypeError("invalid argument for Constraint")

    def __cinit__(self):
        """
        The Cython constructor.

        See :class:`Constraint` for documentation.

        Tests:

            >>> from ppl import Variable
            >>> x = Variable(0)
            >>> x>0   # indirect doctest
            x0>0
        """
        self.thisptr = NULL

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr

    def __hash__(self):
        r"""
        Tests:

        >>> import ppl
        >>> hash(ppl.Variable(0) == 3)
        Traceback (most recent call last):
        ...
        TypeError: Constraint unhashable
        """
        raise TypeError('Constraint unhashable')

    def __repr__(self):
        """
        Return a string representation of the constraint.

        OUTPUT:

        String.

        Examples:

        >>> from ppl import Variable
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
        mpz(3)
        >>> y = Variable(1)
        >>> ineq = 3**50 * y + 2 > 1
        >>> str(ineq.coefficient(y))
        '717897987691852588770249'
        >>> ineq.coefficient(x)
        mpz(0)
        """
        return GMPy_MPZ_From_mpz(self.thisptr.coefficient(v.thisptr[0]).get_mpz_t())

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
        (mpz(3), mpz(5))
        """
        cdef int d = self.space_dimension()
        cdef int i
        coeffs = []
        for i in range(0, d):
            coeffs.append(GMPy_MPZ_From_mpz(self.thisptr.coefficient(PPL_Variable(i)).get_mpz_t()))
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
        mpz(1)
        >>> ineq = 2**66 + y > 0
        >>> str(ineq.inhomogeneous_term())
        '73786976294838206464'
        """
        return GMPy_MPZ_From_mpz(self.thisptr.inhomogeneous_term().get_mpz_t())

    def is_tautological(self):
        r"""
        Test whether ``self`` is a tautological constraint.

        A tautology can have either one of the following forms:

        * an equality: :math:`\sum 0 x_i + 0 = 0`,

        * a non-strict inequality: :math:`\sum 0 x_i + b \geq 0` with `b\geq 0`, or

        * a strict inequality: :math:`\sum 0 x_i + b > 0` with `b> 0`.

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

        * an equality: :math:`\sum 0 x_i + b = 0` with `b\not=0`,

        * a non-strict inequality: :math:`\sum 0 x_i + b \geq 0` with `b< 0`, or

        * a strict inequality: :math:`\sum 0 x_i + b > 0` with `b\leq 0`.

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
        >>> import subprocess, sys
        >>> proc = subprocess.Popen([sys.executable, '-c', cmd], stderr=subprocess.PIPE)
        >>> out, err = proc.communicate()
        >>> print(str(err.decode('ascii')))
        size 4 1 3 2 -1 ...
        <BLANKLINE>
        """
        self.thisptr.ascii_dump()

    def __reduce__(self):
        """
        Pickle object.

        Examples:

        >>> from ppl import Variable
        >>> from pickle import loads, dumps
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

    def permute_space_dimensions(self, cycle):
        """
        Permute the coordinates according to ``cycle``.

        >>> from ppl import Variable
        >>> x = Variable(0); y = Variable(1); z = Variable(2)
        >>> l = 2*x - y + 3*z
        >>> ieq = l >= 5
        >>> ieq.permute_space_dimensions([2, 1, 0])
        >>> ieq
        -x0+3*x1+2*x2-5>=0
        """
        cdef cppvector[PPL_Variable] cpp_cycle
        cdef int i
        for i in cycle:
            cpp_cycle.push_back(PPL_Variable(i))
        self.thisptr.permute_space_dimensions(cpp_cycle)

    def __mod__(self, m):
        r"""
        Return this equality modulo m

        >>> import ppl
        >>> x = ppl.Variable(0)
        >>> y = ppl.Variable(1)
        >>> (2*x + 3*y == 3) % 5
        2*x0+3*x1-3==0 (mod 5)

        >>> (x <= 3) % 5
        Traceback (most recent call last):
        ...
        ValueError: PPL::Congruence::Congruence(c, r):
        constraint c must be an equality.
        """
        return Congruence(self, m)


def inequality(expression):
    """
    Construct an inequality.

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


def strict_inequality(expression):
    """
    Construct a strict inequality.

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


def equation(expression):
    """
    Construct an equation.

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
        >>> cs.insert( 6*x < 3*y )
        >>> cs.insert( x >= 2*x-7*y )
        >>> cs
        Constraint_System {5*x0-2*x1>0, ...}
        >>> cs[0]
        5*x0-2*x1>0
    """
    def __init__(self, arg=None):
        if arg is None:
            self.thisptr = new PPL_Constraint_System()
        elif isinstance(arg, Constraint):
            g = <Constraint>arg
            self.thisptr = new PPL_Constraint_System(g.thisptr[0])
        elif isinstance(arg, Constraint_System):
            gs = <Constraint_System>arg
            self.thisptr = new PPL_Constraint_System(gs.thisptr[0])
        elif isinstance(arg, (list, tuple)):
            self.thisptr = new PPL_Constraint_System()
            for constraint in arg:
                self.insert(constraint)
        else:
            raise TypeError('cannot initialize from {!r}'.format(arg))

    def __cinit__(self, arg=None):
        """
        The Cython constructor.

        See :class:`Constraint_System` for documentation.

        Tests:

        >>> from ppl import Constraint_System
        >>> Constraint_System()
        Constraint_System {}
        """
        self.thisptr = NULL

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr

    def __hash__(self):
        r"""
        Tests:

        >>> import ppl
        >>> hash(ppl.Constraint_System())
        Traceback (most recent call last):
        ...
        TypeError: Constraint_System unhashable
        """
        raise TypeError('Constraint_System unhashable')

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
        self.thisptr.insert(c.thisptr[0])

    def empty(self):
        """
        Return ``True`` if and only if ``self`` has no constraints.

        OUTPUT:

        Boolean.

        Examples:

        >>> from ppl import Variable, Constraint_System
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
        >>> import subprocess, sys
        >>> proc = subprocess.Popen([sys.executable, '-c', cmd], stderr=subprocess.PIPE)
        >>> out, err = proc.communicate()
        >>> print(str(err.decode('ascii')))
        topology NOT_NECESSARILY_CLOSED
        ...
        <BLANKLINE>
        """
        self.thisptr.ascii_dump()

    def __len__(self):
        """
        Return the number of constraints in the system.

        Examples:

        >>> from ppl import Variable, Constraint_System
        >>> x = Variable(0)
        >>> cs = Constraint_System( [x>0, x<1] )
        >>> len(cs)
        2
        >>> cs.insert(2*x < 3)
        >>> len(cs)
        3
        >>> cs.clear()
        >>> len(cs)
        0
        """
        cdef Py_ssize_t length = 0
        cdef PPL_Constraint_System_iterator *csi_ptr = new PPL_Constraint_System_iterator(self.thisptr[0].begin())
        while csi_ptr[0] != self.thisptr[0].end():
            length += 1
            csi_ptr[0].inc(1)
        del csi_ptr
        return length

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
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System( 5*x < 2*y )
        >>> cs.insert( 6*x-y == 0 )
        >>> cs.insert( x >= 2*x-7*y )
        >>> next(cs.__iter__())
        -5*x0+2*x1>0
        >>> list(cs)
        [-5*x0+2*x1>0, 6*x0-x1==0, -x0+7*x1>=0]
        >>> x = Variable(0)
        >>> cs = Constraint_System( x > 0 )
        >>> cs.insert ( 2*x <= -3)
        >>> it = cs.__iter__()
        >>> next(it)
        x0>0
        >>> next(it)
        -2*x0-3>=0
        >>> next(it)
        Traceback (most recent call last):
        ...
        StopIteration
        """
        cdef PPL_Constraint_System_iterator *csi_ptr = new PPL_Constraint_System_iterator(self.thisptr[0].begin())
        try:
            while csi_ptr[0] != self.thisptr[0].end():
                yield _wrap_Constraint(deref(csi_ptr[0].inc(1)))
        finally:
            del csi_ptr

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
        s += ', '.join([repr(c) for c in self])
        s += '}'
        return s

    def __reduce__(self):
        """
        Pickle object.

        Tests:

        >>> from ppl import Constraint_System, Variable
        >>> from pickle import loads, dumps
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Constraint_System([3*x+2*y+1 < 3, 0*x>x+1]);  cs
        Constraint_System {-3*x0-2*x1+2>0, -x0-1>0}
        >>> loads(dumps(cs))
        Constraint_System {-3*x0-2*x1+2>0, -x0-1>0}
        """
        return (Constraint_System, (tuple(self),))


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
    ...       s = ""
    ...       for j, rel_j in enumerate(rels):
    ...           s=s+str(int(rel_i.implies(rel_j))) + ' '
    ...       print(" ".join(str(int(rel_i.implies(rel_j))) for j, rel_j in enumerate(rels)))
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

    def __hash__(self):
        r"""
        Tests:

        >>> import ppl
        >>> hash(ppl.Poly_Con_Relation.saturates())
        Traceback (most recent call last):
        ...
        TypeError: Poly_Con_Relation unhashable
        """
        raise TypeError('Poly_Con_Relation unhashable')

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
        >>> import subprocess, sys
        >>> proc = subprocess.Popen([sys.executable, '-c', cmd], stderr=subprocess.PIPE)
        >>> out, err = proc.communicate()
        >>> print(str(err.decode('ascii')))
        NOTHING
        """
        self.thisptr.ascii_dump()

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


cdef _make_Constraint_from_richcmp(lhs_, rhs_, op):
    cdef Linear_Expression lhs = Linear_Expression(lhs_)
    cdef Linear_Expression rhs = Linear_Expression(rhs_)
    if op == Py_LT:
        return _wrap_Constraint(lhs.thisptr[0] < rhs.thisptr[0])
    elif op == Py_LE:
        return _wrap_Constraint(lhs.thisptr[0] <= rhs.thisptr[0])
    elif op == Py_EQ:
        return _wrap_Constraint(lhs.thisptr[0] == rhs.thisptr[0])
    elif op == Py_GT:
        return _wrap_Constraint(lhs.thisptr[0] > rhs.thisptr[0])
    elif op == Py_GE:
        return _wrap_Constraint(lhs.thisptr[0] >= rhs.thisptr[0])
    elif op == Py_NE:
        raise NotImplementedError
    else:
        assert(False)


cdef _wrap_Constraint(PPL_Constraint constraint):
    """
    Wrap a C++ ``PPL_Constraint`` into a Cython ``Constraint``.
    """
    cdef Constraint c = Constraint.__new__(Constraint)
    c.thisptr = new PPL_Constraint(constraint)
    return c


cdef _wrap_Constraint_System(PPL_Constraint_System constraint_system):
    """
    Wrap a C++ ``PPL_Constraint_System`` into a Cython ``Constraint_System``.
    """
    cdef Constraint_System cs = Constraint_System.__new__(Constraint_System)
    cs.thisptr = new PPL_Constraint_System(constraint_system)
    return cs


cdef _wrap_Poly_Con_Relation(PPL_Poly_Con_Relation relation):
    """
    Wrap a C++ ``PPL_Poly_Con_Relation`` into a Cython ``Poly_Con_Relation``.
    """
    cdef Poly_Con_Relation rel = Poly_Con_Relation(True)
    rel.thisptr = new PPL_Poly_Con_Relation(relation)
    return rel
