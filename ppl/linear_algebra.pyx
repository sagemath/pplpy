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

cimport cython

from gmpy2 cimport import_gmpy2, mpz, GMPy_MPZ_From_mpz, MPZ_Check
from .constraint cimport _make_Constraint_from_richcmp

# initialize gmpy2 C API
import_gmpy2()

####################################################
# PPL can use floating-point arithmetic to compute integers
cdef extern from "ppl.hh" namespace "Parma_Polyhedra_Library":
    cdef void set_rounding_for_PPL()
    cdef void restore_pre_PPL_rounding()

# but with PPL's rounding the gsl will be very unhappy; must turn off!
restore_pre_PPL_rounding()

cdef PPL_Coefficient PPL_Coefficient_from_pyobject(c) except *:
    cdef mpz coeff

    if MPZ_Check(c):
        coeff = <mpz> c
    elif isinstance(c, (int, str)):
        coeff = mpz(c)
    else:
        coeff = mpz(c)
        if coeff != c or c != coeff:
            raise TypeError('ppl coefficients must be integral')

    return PPL_Coefficient(coeff.z)


@cython.freelist(128)
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

    def __hash__(self):
        r"""
        Tests:

        >>> import ppl
        >>> hash(ppl.Variable(12))
        Traceback (most recent call last):
        ...
        TypeError: Variable unhashable
        """
        raise TypeError('Variable unhashable')

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

    def space_dimension(self):
        r"""
        Return the dimension of the vector space enclosing ``self``.

        OUTPUT:

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

        >>> from gmpy2 import mpz
        >>> x + mpz(3)
        x0+3
        >>> mpz(-5) + y
        x1-5

        >>> x + 1.5
        Traceback (most recent call last):
        ...
        TypeError: ppl coefficients must be integral
        >>> 1.5 + x
        Traceback (most recent call last):
        ...
        TypeError: ppl coefficients must be integral
        """
        return Linear_Expression(self) + Linear_Expression(other)

    def __radd__(self, other):
        return Linear_Expression(self) + Linear_Expression(other)

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
        return Linear_Expression(self) - Linear_Expression(other)

    def __rsub__(self, other):
        return Linear_Expression(other) - Linear_Expression(self)

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

        >>> 1.5 * x
        Traceback (most recent call last):
        ...
        TypeError: ppl coefficients must be integral
        >>> x * 1.5
        Traceback (most recent call last):
        ...
        TypeError: ppl coefficients must be integral
        """
        if isinstance(self, Variable):
            return Linear_Expression(self) * other
        else:
            # NOTE: this code path will only be executed when compiled with cython < 3.0.0
            return Linear_Expression(other) * self

    def __rmul__(self, other):
        return Linear_Expression(self) * other

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
### Variables_Set ##################################
####################################################
cdef class Variables_Set(object):
    r"""
    Wrapper for PPL's ``Variables_Set`` class.

    A set of variables' indexes.

    Examples:

    Build the empty set of variable indexes:

    >>> from ppl import Variable, Variables_Set
    >>> Variables_Set()
    Variables_Set of cardinality 0

    Build the singleton set of indexes containing the index of the variable:

    >>> v123 = Variable(123)
    >>> Variables_Set(v123)
    Variables_Set of cardinality 1

    Build the set of variables' indexes in the range from one variable to
    another variable:

    >>> v127 = Variable(127)
    >>> Variables_Set(v123,v127)
    Variables_Set of cardinality 5

    You can alternatively use integers:

    >>> Variables_Set(23)
    Variables_Set of cardinality 1
    >>> Variables_Set(10, 14)
    Variables_Set of cardinality 5
    """
    def __cinit__(self, *args):
        """
        The Cython constructor.

        See :class:`Variables_Set` for documentation.

        Tests:

        >>> from ppl import Variable, Variables_Set
        >>> Variables_Set()
        Variables_Set of cardinality 0
        """
        cdef Variable arg0, arg1
        if len(args) == 0:
            self.thisptr = new PPL_Variables_Set()
        elif len(args) == 1:
            if type(args[0]) is Variable:
                arg0 = <Variable> args[0]
            else:
                arg0 = Variable(args[0])
            self.thisptr = new PPL_Variables_Set(arg0.thisptr[0])
        elif len(args) == 2:
            if type(args[0]) is Variable:
                arg0 = <Variable> args[0]
            else:
                arg0 = Variable(args[0])

            if type(args[1]) is Variable:
                arg1 = <Variable> args[1]
            else:
                arg1 = Variable(args[1])

            self.thisptr = new PPL_Variables_Set(arg0.thisptr[0], arg1.thisptr[0])

    def __hash__(self):
        r"""
        Tests:

        >>> import ppl
        >>> hash(ppl.Variables_Set())
        Traceback (most recent call last):
        ...
        TypeError: Variables_Set unhashable
        """
        raise TypeError('Variables_Set unhashable')

    def __dealloc__(self):
        """
        The Cython destructor
        """
        del self.thisptr

    def space_dimension(self):
        r"""
        Returns the dimension of the smallest vector space enclosing all the variables whose indexes are in the set.

        OUTPUT:

        Integer.

        EXAMPLES::

            >>> from ppl import Variable, Variables_Set
            >>> v123 = Variable(123)
            >>> S = Variables_Set(v123)
            >>> S.space_dimension()
            124
        """
        return self.thisptr.space_dimension()

    def insert(self, Variable v):
        r"""
        Inserts the index of variable `v` into the set.

        EXAMPLES::

            >>> from ppl import Variable, Variables_Set
            >>> S = Variables_Set()
            >>> v123 = Variable(123)
            >>> S.insert(v123)
            >>> S.space_dimension()
            124
        """
        self.thisptr.insert(v.thisptr[0])

    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        TODO: rewrite examples
        EXAMPLES::

            >>> cmd  = 'from ppl import Variable, Variables_Set\n'
            >>> cmd += 'v123 = Variable(123)\n'
            >>> cmd += 'S = Variables_Set(v123)\n'
            >>> cmd += 'S.ascii_dump()\n'
            >>> import subprocess
            >>> import sys
            >>> proc = subprocess.Popen([sys.executable, '-c', cmd], stderr=subprocess.PIPE)
            >>> out, err = proc.communicate()
            >>> print(str(err.decode('ascii')))
            <BLANKLINE>
            variables( 1 )
            123
        """
        self.thisptr.ascii_dump()

    def __repr__(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            >>> from ppl import Variable, Variables_Set
            >>> S = Variables_Set()
            >>> S.__repr__()
            'Variables_Set of cardinality 0'
        """
        return 'Variables_Set of cardinality {}'.format(self.thisptr.size())

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

    - ``a`` -- either a dictionary whose indices are space dimension and
      values are coefficients or an iterable coefficients (e.g. a list or
      tuple).

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

    >>> Linear_Expression({1: -3, 7: 1}, 0)
    -3*x1+x7
    >>> Linear_Expression([1, 2, 3, 4], 5)
    x0+2*x1+3*x2+4*x3+5
    >>> Linear_Expression(10)
    10
    >>> Linear_Expression()
    0
    >>> Linear_Expression({}, 2)
    2
    >>> Linear_Expression([], 2)
    2

    >>> Linear_Expression(10).inhomogeneous_term()
    mpz(10)
    >>> x = Variable(123)
    >>> expr = x+1
    >>> expr
    x123+1
    >>> expr.coefficient(x)
    mpz(1)
    >>> expr.coefficient(Variable(124))
    mpz(0)

    >>> from gmpy2 import mpz, mpq
    >>> Linear_Expression(mpz(3))
    3
    >>> Linear_Expression([mpz(5), mpz(2)], mpz(-2))
    5*x0+2*x1-2

    String, rationals and floating point types are accepted as long as they
    represent exact integers:

    >>> Linear_Expression(('4', 1), 2)
    4*x0+x1+2
    >>> Linear_Expression((4, 1.0, mpq('4/2')), 2.0)
    4*x0+x1+2*x2+2
    >>> Linear_Expression(('1.5',), 0)
    Traceback (most recent call last):
    ...
    ValueError: invalid digits
    >>> Linear_Expression((mpq('3/2'),), 0)
    Traceback (most recent call last):
    ...
    TypeError: ppl coefficients must be integral
    >>> Linear_Expression((1, 2.1, 1), 1)
    Traceback (most recent call last):
    ...
    TypeError: ppl coefficients must be integral
    >>> Linear_Expression(mpq('1/2'))
    Traceback (most recent call last):
    ...
    TypeError: ppl coefficients must be integral
    >>> Linear_Expression('I am a linear expression')
    Traceback (most recent call last):
    ...
    ValueError: invalid digits
    >>> Linear_Expression(('I','am','a','linear','expression'))
    Traceback (most recent call last):
    ...
    TypeError: mpz() requires numeric or string argument
    """
    def __init__(self, *args):
        """
        The Cython constructor.

        See :class:`Linear_Expression` for documentation.

        Tests:

        >>> from ppl import Linear_Expression
        >>> Linear_Expression(10)   # indirect doctest
        10
        """
        cdef size_t i
        if len(args) == 2:
            a = args[0]
            b = args[1]
            self.thisptr = new PPL_Linear_Expression()
            if isinstance(a, dict):
                if a:
                    self.thisptr.set_space_dimension(1 + max(a))
                    for i, coeff in a.items():
                        self.thisptr.set_coefficient(PPL_Variable(i), PPL_Coefficient_from_pyobject(coeff))
            else:
                self.thisptr.set_space_dimension(len(a))
                for i, coeff in enumerate(a):
                    self.thisptr.set_coefficient(PPL_Variable(i), PPL_Coefficient_from_pyobject(coeff))
            self.thisptr.set_inhomogeneous_term(PPL_Coefficient_from_pyobject(b))
            return
        elif len(args) == 1:
            arg = args[0]
            if isinstance(arg, Variable):
                v = <Variable>arg
                self.thisptr = new PPL_Linear_Expression(v.thisptr[0])
                return
            if isinstance(arg, Linear_Expression):
                e = <Linear_Expression>arg
                self.thisptr = new PPL_Linear_Expression(e.thisptr[0])
                return
            self.thisptr = new PPL_Linear_Expression(PPL_Coefficient_from_pyobject(arg))
        elif len(args) == 0:
            self.thisptr = new PPL_Linear_Expression()
            return
        else:
            raise ValueError("Cannot initialize with more than 2 arguments.")

    def __dealloc__(self):
        """
        The Cython destructor.
        """
        del self.thisptr

    def __hash__(self):
        r"""
        Tests:

        >>> import ppl
        >>> hash(ppl.Linear_Expression(10))
        Traceback (most recent call last):
        ...
        TypeError: Linear_Expression unhashable
        """
        raise TypeError('Linear_Expression unhashable')

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

    def set_space_dimension(self, PPL_dimension_type n):
        r"""
        Set the dimension of the ambient space to ``n``

        Examples:

        >>> import ppl
        >>> L = ppl.Variable(0) + ppl.Variable(3)
        >>> L.space_dimension()
        4
        >>> L.set_space_dimension(6)
        >>> L.space_dimension()
        6

        >>> L = ppl.Variable(5) - ppl.Variable(2)
        >>> L.set_space_dimension(3)
        >>> L
        -x2
        """
        self.thisptr.set_space_dimension(n)

    def swap_space_dimensions(self, v1, v2):
        r"""
        Swaps the coefficients of ``v1`` and ``v2``.

        INPUT:

        - ``v1``, ``v2`` - variables or indices of variables

        Examples:

        >>> import ppl
        >>> L = ppl.Variable(1) - 3 * ppl.Variable(3)
        >>> L.swap_space_dimensions(ppl.Variable(1), ppl.Variable(3))
        >>> L
        -3*x1+x3

        >>> L = ppl.Variable(1) - 3 * ppl.Variable(3)
        >>> L.swap_space_dimensions(1, 3)
        >>> L
        -3*x1+x3
        """
        cdef Variable vv1, vv2
        if type(v1) is Variable:
            vv1 = <Variable> v1
        else:
            vv1 = Variable(v1)
        if type(v2) is Variable:
            vv2 = <Variable> v2
        else:
            vv2 = Variable(v2)
        self.thisptr.swap_space_dimensions(vv1.thisptr[0], vv2.thisptr[0])

    def shift_space_dimensions(self, v, PPL_dimension_type n):
        r"""
        Shift by ``n`` the coefficients of variables starting from the
        coefficient of ``v``.

        This increases the space dimension by ``n``.

        Examples:

        >>> import ppl
        >>> L = ppl.Variable(0) + 13 * ppl.Variable(2) + 5 * ppl.Variable(7)
        >>> L
        x0+13*x2+5*x7
        >>> L.shift_space_dimensions(ppl.Variable(2), 2)
        >>> L
        x0+13*x4+5*x9
        >>> L.shift_space_dimensions(ppl.Variable(7), 3)
        >>> L
        x0+13*x4+5*x12
        """
        cdef Variable vv
        if type(v) is Variable:
            vv = <Variable> v
        else:
            vv = Variable(v)
        self.thisptr.shift_space_dimensions(vv.thisptr[0], n)

    def remove_space_dimensions(self, Variables_Set V):
        r"""
        Removes the dimension specified by the set of variables ``V``.

        See :class:`Variables_Set` to construct set of variables.

        Examples:

        >>> import ppl
        >>> L = sum(i * ppl.Variable(i) for i in range(10))
        >>> L
        x1+2*x2+3*x3+4*x4+5*x5+6*x6+7*x7+8*x8+9*x9
        >>> L.remove_space_dimensions(ppl.Variables_Set(3,5))
        >>> L
        x1+2*x2+6*x3+7*x4+8*x5+9*x6
        """
        self.thisptr.remove_space_dimensions(V.thisptr[0])

    def coefficient(self, v):
        """
        Return the coefficient of the variable ``v``.

        INPUT:

        - ``v`` -- a :class:`Variable`.

        OUTPUT:

        A gmpy2 integer.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0)
        >>> e = 3*x+1
        >>> e.coefficient(x)
        mpz(3)
        >>> e.coefficient(Variable(13))
        mpz(0)
        """
        cdef Variable vv
        if type(v) is Variable:
            vv = <Variable> v
        else:
            vv = Variable(v)
        return GMPy_MPZ_From_mpz(self.thisptr.coefficient(vv.thisptr[0]).get_mpz_t())

    def coefficients(self):
        """
        Return the coefficients of the linear expression.

        OUTPUT:

        A tuple of gmpy2 integers of length :meth:`space_dimension`.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0);  y = Variable(1)
        >>> e = 3*x+5*y+1
        >>> e.coefficients()
        (mpz(3), mpz(5))
        """
        cdef int d = self.space_dimension()
        cdef int i
        cdef list coeffs = [None]*d
        for i in range(d):
            coeffs[i] = GMPy_MPZ_From_mpz(self.thisptr.coefficient(PPL_Variable(i)).get_mpz_t())
        return tuple(coeffs)

    def set_coefficient(self, i, v):
        """
        Set the ``i``-th coefficient to ``v``.

        INPUT:

        - ``i`` - variable or variable index

        - ``v`` - integer

        Examples:

        >>> from ppl import Variable
        >>> L = Variable(0) + 3 * Variable(1)
        >>> L.set_coefficient(1, -5)
        >>> L.set_coefficient(7, 3)
        >>> L
        x0-5*x1
        """
        cdef Variable ii
        if type(i) is Variable:
            ii = <Variable> i
        else:
            ii = Variable(i)

        cdef PPL_Coefficient vv = PPL_Coefficient_from_pyobject(v)

        self.thisptr.set_coefficient(ii.thisptr[0], vv)

    def inhomogeneous_term(self):
        """
        Return the inhomogeneous term of the linear expression.

        OUTPUT:

        Integer.

        Examples:

        >>> from ppl import Linear_Expression
        >>> Linear_Expression(10).inhomogeneous_term()
        mpz(10)
        """
        return GMPy_MPZ_From_mpz(self.thisptr.inhomogeneous_term().get_mpz_t())

    def set_inhomogeneous_term(self, v):
        """
        Set the inhomogeneous term of this linear expression.

        Examples:

        >>> from ppl import Linear_Expression
        >>> L = Linear_Expression()
        >>> L.set_inhomogeneous_term(-1313958534578713747)
        >>> L.inhomogeneous_term()
        mpz(-1313958534578713747)
        """
        cdef PPL_Coefficient vv = PPL_Coefficient_from_pyobject(v)
        self.thisptr.set_inhomogeneous_term(vv)

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

    def is_equal_to(self, Linear_Expression other):
        """
        Test equality with another linear expression.

        OUTPUT: boolean

        Examples:

        >>> from ppl import Variable
        >>> L1 = Variable(0) + 2 * Variable(3)
        >>> L2 = Variable(0) + 2 * Variable(3)
        >>> L3 = Variable(0) - Variable(2)
        >>> L1.is_equal_to(L2)
        True
        >>> L1.is_equal_to(L3)
        False
        """
        return self.thisptr.is_equal_to(other.thisptr[0])

    def ascii_dump(self):
        r"""
        Write an ASCII dump to stderr.

        Examples:

        >>> cmd  = 'from ppl import Linear_Expression, Variable\n'
        >>> cmd += 'x = Variable(0)\n'
        >>> cmd += 'y = Variable(1)\n'
        >>> cmd += 'e = 3*x+2*y+1\n'
        >>> cmd += 'e.ascii_dump()\n'
        >>> from subprocess import Popen, PIPE
        >>> import sys
        >>> proc = Popen([sys.executable, '-c', cmd], stdout=PIPE, stderr=PIPE)
        >>> out, err = proc.communicate()
        >>> len(out) == 0
        True
        >>> len(err) > 0
        True
        """
        self.thisptr.ascii_dump()

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
        for i in range(self.space_dimension()):
            x = Variable(i)
            coeff = self.coefficient(x)
            if coeff == 0:
                continue
            if first and coeff == 1:
                s += '%r' % x
                first = False
            elif first and coeff == -1:
                s += '-%r' % x
                first = False
            elif first and coeff != 1:
                s += '%d*%r' % (coeff, x)
                first = False
            elif coeff == 1:
                s += '+%r' % x
            elif coeff == -1:
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
        >>> 9 + x + y + (1 + x) + y + y
        2*x0+3*x1+10

        >>> from gmpy2 import mpz
        >>> mpz(3) + x + mpz(5) + y + mpz(7)
        x0+x1+15
        >>> from gmpy2 import mpz
        >>> x + mpz(5)
        x0+5
        """
        cdef Linear_Expression lhs, rhs

        if isinstance(self, Linear_Expression):
            lhs = <Linear_Expression> self
        else:
            # NOTE: this code path will only be executed when compiled with cython < 3.0.0
            lhs = Linear_Expression(self)

        if isinstance(other, Linear_Expression):
            rhs = <Linear_Expression> other
        else:
            rhs = Linear_Expression(other)

        cdef Linear_Expression result = Linear_Expression()
        result.thisptr[0] = lhs.thisptr[0] + rhs.thisptr[0]
        return result

    def __radd__(self, other):
        r"""
        Add ``self`` and ``other``.
        """
        cdef Linear_Expression lhs, rhs

        lhs = <Linear_Expression> self

        if isinstance(other, Linear_Expression):
            rhs = <Linear_Expression> other
        else:
            rhs = Linear_Expression(other)

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

        >>> from gmpy2 import mpz
        >>> mpz(5)-x-(mpz(3)-y)-x-mpz(7)
        -2*x0+x1-5
        """
        cdef Linear_Expression lhs, rhs

        if isinstance(self, Linear_Expression):
            lhs = <Linear_Expression> self
        else:
            lhs = Linear_Expression(self)

        if isinstance(other, Linear_Expression):
            rhs = <Linear_Expression> other
        else:
            rhs = Linear_Expression(other)

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
        >>> 8 * (x + 1)
        8*x0+8
        >>> y * 8
        8*x1
        >>> 2**128 * x
        340282366920938463463374607431768211456*x0

        >>> from gmpy2 import mpz
        >>> mpz(3) * x * mpz(5)
        15*x0
        """
        cdef Linear_Expression e
        cdef c

        if isinstance(self, Linear_Expression):
            e = <Linear_Expression> self
            c = other
        else:
            # NOTE: this code path will only be executed when compiled with cython < 3.0.0
            e = <Linear_Expression> other
            c = self

        cdef PPL_Coefficient cc = PPL_Coefficient_from_pyobject(c)
        cdef Linear_Expression result = Linear_Expression()
        result.thisptr[0] = e.thisptr[0] * cc
        return result

    def __rmul__(self, other):
        r"""
        Multiply ``self`` with ``other``
        """
        cdef Linear_Expression e
        cdef c
        e = <Linear_Expression> self
        c = other

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
        >>> from pickle import loads, dumps
        >>> le = loads(dumps(Linear_Expression([1,2,3],4)))
        >>> le.coefficients() == (1,2,3)
        True
        >>> le.inhomogeneous_term() == 4
        True
        """
        return (Linear_Expression, (self.coefficients(), self.inhomogeneous_term()))

    def permute_space_dimensions(self, cycle):
        r"""
        Permute the coordinates according to ``cycle``.

        Examples:

        >>> from ppl import Variable
        >>> x = Variable(0); y = Variable(1); z = Variable(2)
        >>> l = 2*x - y + 3*z
        >>> l.permute_space_dimensions([0, 2])
        >>> l
        3*x0-x1+2*x2
        >>> l.permute_space_dimensions([1, 0, 2])
        >>> l
        -x0+2*x1+3*x2
        """
        cdef cppvector[PPL_Variable] cpp_cycle
        cdef PPL_dimension_type i
        for i in cycle:
            cpp_cycle.push_back(PPL_Variable(i))
        self.thisptr.permute_space_dimensions(cpp_cycle)
