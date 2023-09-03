# distutils: language = c++
# distutils: libraries = gmp gmpxx ppl m
#*****************************************************************************
#       Copyright (C) 2020 Vincent Delecroix <vincent.delecroix@labri.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cython.operator cimport dereference as deref

from gmpy2 cimport GMPy_MPZ_From_mpz, mpz, import_gmpy2

from .constraint cimport Constraint
from .linear_algebra cimport Variable, Linear_Expression, PPL_Coefficient_from_pyobject

import_gmpy2()


def _dummy():
    raise ValueError


cdef class Congruence(object):
    r"""
    Wrapper for PPL's ``Congruence`` class.

    >>> import ppl
    >>> x = ppl.Variable(0)
    >>> y = ppl.Variable(1)
    >>> ppl.Congruence(x + 2*y - 1, 7)
    x0+2*x1-1==0 (mod 7)
    >>> ppl.Congruence(x + y == 2, 5)
    x0+x1-2==0 (mod 5)
    """
    def __init__(self, arg=None, mod=None):
        if arg is None:
            self.thisptr = new PPL_Congruence()
        elif isinstance(arg, Congruence):
            self.thisptr = new PPL_Congruence((<Congruence> arg).thisptr[0])
        elif isinstance(arg, Variable):
            self.thisptr = new PPL_Congruence()
            if mod is None:
                mod = mpz(1)
            self.thisptr[0] = modulo(PPL_Linear_Expression((<Variable> arg).thisptr[0]),
                                     PPL_Coefficient_from_pyobject(mod))
        elif isinstance(arg, Linear_Expression):
            self.thisptr = new PPL_Congruence()
            if mod is None:
                mod = mpz(1)
            self.thisptr[0] = modulo((<Linear_Expression> arg).thisptr[0],
                                     PPL_Coefficient_from_pyobject(mod))
        elif isinstance(arg, Constraint):
            self.thisptr = new PPL_Congruence((<Constraint> arg).thisptr[0])
            if mod is None:
                mod = mpz(1)
            self.thisptr.set_modulus(PPL_Coefficient_from_pyobject(mod))
        else:
            raise TypeError("invalid argument for Congruence")

    def __cinit__(self):
        self.thisptr = NULL

    def __dealloc__(self):
        del self.thisptr

    def is_equal(self, Congruence other):
        r"""
        Return whether ``self`` is equal to ``other``.

        Examples:

        >>> import ppl
        >>> x = ppl.Variable(0)
        >>> Congruence(x + 1, 3).is_equal(Congruence(x + 1, 3))
        True
        >>> Congruence(x, 3).is_equal(Congruence(x - 2, 3))
        False
        >>> Congruence(x, 3).is_equal(Congruence(x , 2))
        False
        """
        return (<Congruence> self).thisptr[0] == (<Congruence> other).thisptr[0]

    def coefficient(self, v):
        r"""
        Return the coefficient ``v`` of this congruence.

        Examples:

        >>> import ppl
        >>> x = ppl.Variable(0)
        >>> y = ppl.Variable(1)
        >>> c = (2 * x + y == 1) % 12
        >>> c.coefficient(0)
        mpz(2)
        >>> c.coefficient(x)
        mpz(2)

        Note that contrarily to :class:`Linear_Expression` the congruences raise
        an error when trying to access a coefficient beyond the dimension

        >>> c.coefficient(ppl.Variable(13))
        Traceback (most recent call last):
        ...
        ValueError: PPL::Congruence::coefficient(v):
        this->space_dimension() == 2, v.space_dimension() == 14.
        """
        cdef Variable vv
        if type(v) is Variable:
            vv = <Variable> v
        else:
            vv = Variable(v)
        return GMPy_MPZ_From_mpz(self.thisptr.coefficient(vv.thisptr[0]).get_mpz_t())

    def coefficients(self):
        r"""
        Return the coefficients of this congruence as a tuple.

        Examples:

        >>> import ppl
        >>> x = ppl.Variable(0)
        >>> y = ppl.Variable(1)
        >>> t = ppl.Variable(3)
        >>> c = ppl.Congruence(x + 2*y + t - 1, 7)
        >>> c.coefficients()
        (mpz(1), mpz(2), mpz(0), mpz(1))
        """
        cdef int d = self.space_dimension()
        cdef int i
        coeffs = []
        for i in range(0, d):
            coeffs.append(GMPy_MPZ_From_mpz(self.thisptr.coefficient(PPL_Variable(i)).get_mpz_t()))
        return tuple(coeffs)

    def modulus(self):
        r"""
        Return the modulus of this congruence.

        Examples:

        >>> import ppl
        >>> x = ppl.Variable(0)
        >>> y = ppl.Variable(1)
        >>> c = ppl.Congruence(x + 2*y - 3, 7)
        >>> c.modulus()
        mpz(7)
        """
        return GMPy_MPZ_From_mpz(self.thisptr.modulus().get_mpz_t())

    def inhomogeneous_term(self):
        r"""
        Return the inhomogeneous term of this congruence.

        Examples:

        >>> import ppl
        >>> x = ppl.Variable(0)
        >>> y = ppl.Variable(1)
        >>> c = ppl.Congruence(x + 2*y - 3, 7)
        >>> c.inhomogeneous_term()
        mpz(-3)
        """
        return GMPy_MPZ_From_mpz(self.thisptr.inhomogeneous_term().get_mpz_t())

    def space_dimension(self):
        return self.thisptr.space_dimension()

    def __repr__(self):
        s = ''
        v = [Variable(i) for i in range(self.space_dimension())]
        e = sum(self.coefficient(x)*x for x in v)
        e += self.inhomogeneous_term()
        s = repr(e) + '==0 (mod %s)' % self.modulus()
        return s

    def __reduce__(self):
        r"""
        >>> from ppl import Variable, Congruence, Congruence_System
        >>> from pickle import loads, dumps
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> z = Variable(2)
        >>> c1 = Congruence(2*x + 3*y - 5*z == 3, 12)
        >>> c1
        2*x0+3*x1-5*x2-3==0 (mod 12)
        >>> loads(dumps(c1))
        2*x0+3*x1-5*x2-3==0 (mod 12)
        >>> assert c1.is_equal(loads(dumps(c1)))
        """
        le = Linear_Expression(self.coefficients(), self.inhomogeneous_term())
        return (congruence, (le == 0, self.modulus()))

cdef class Congruence_System(object):
    r"""
    Wrapper for PPL's ``Congruence_System`` class.

    >>> from ppl import Variable, Congruence, Congruence_System
    >>> x = Variable(0)
    >>> y = Variable(1)
    >>> z = Variable(2)
    >>> c1 = Congruence(2*x + 3*y - 5*z == 0, 12)
    >>> c2 = Congruence(4*x + y == 5, 18)
    >>> C = Congruence_System()
    >>> C.insert(c1)
    >>> C.insert(c2)
    >>> C
    Congruence_System {2*x0+3*x1-5*x2==0 (mod 12), 4*x0+x1+13==0 (mod 18)}
    """
    def __init__(self, arg=None):
        if arg is None:
            self.thisptr = new PPL_Congruence_System()
        elif isinstance(arg, Congruence):
            c = <Congruence> arg
            self.thisptr = new PPL_Congruence_System(c.thisptr[0])
        elif isinstance(arg, Congruence_System):
            cs = <Congruence_System> arg
            self.thisptr = new PPL_Congruence_System(cs.thisptr[0])
        elif isinstance(arg, (list, tuple)):
            self.thisptr = new PPL_Congruence_System()
            for congruence in arg:
                self.insert(congruence)

    def ascii_dump(self):
        r"""
        Write an ASCII dump of this congruence to stderr.

        Examples:

        >>> cmd  = 'from ppl import Variable, Congruence_System\n'
        >>> cmd += 'x = Variable(0)\n'
        >>> cmd += 'y = Variable(1)\n'
        >>> cmd += 'cs = Congruence_System( (3*x == 2*y+1) % 7 )\n'
        >>> cmd += 'cs.ascii_dump()\n'
        >>> import subprocess, sys
        >>> proc = subprocess.Popen([sys.executable, '-c', cmd], stderr=subprocess.PIPE)
        >>> out, err = proc.communicate()
        >>> print(str(err.decode('ascii')))
        1 x 2 SPARSE
        size 3 6 3 -2 m 7
        <BLANKLINE>
        """
        self.thisptr.ascii_dump()

    def __len__(self):
        r"""
        Return the number of congruences in the system.

        Examples:

        >>> from ppl import Variable, Congruence_System
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Congruence_System( [(x == 3) % 15, (2*y == 7) % 12])
        >>> len(cs)
        2
        >>> cs.insert((2*x + y == 12) % 22)
        >>> len(cs)
        3
        >>> cs.clear()
        >>> len(cs)
        0
        """
        cdef Py_ssize_t l = 0
        cdef PPL_Congruence_System_iterator * csi_ptr = new PPL_Congruence_System_iterator(self.thisptr[0].begin())

        while csi_ptr[0] != self.thisptr[0].end():
            l += 1
            csi_ptr[0].inc(1)
        del csi_ptr
        return l

    def __iter__(self):
        cdef PPL_Congruence_System_iterator * csi_ptr = new PPL_Congruence_System_iterator(self.thisptr[0].begin())

        try:
            while csi_ptr[0] != self.thisptr[0].end():
                yield _wrap_Congruence(deref(csi_ptr[0].inc(1)))
        finally:
            del csi_ptr

    def __repr__(self):
        s = 'Congruence_System {'
        s += ', '.join([repr(c) for c in self])
        s += '}'
        return s

    def insert(self, Congruence c):
        r"""
        Insert the congruence ``c`` into the congruence system.

        Examples:

        >>> from ppl import Variable, Congruence_System
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Congruence_System()
        >>> cs.insert((3 * x + 2 == 0) % 13)
        >>> cs
        Congruence_System {3*x0+2==0 (mod 13)}
        """
        self.thisptr.insert(c.thisptr[0])

    def clear(self):
        r"""
        Remove all congruences from this congruence system.

        >>> from ppl import Variable, Congruence_System
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Congruence_System([(3*x+2*y+1 == 3) % 15, (x+2*y + 7 == 0) % 22])
        >>> len(cs)
        2
        >>> cs.clear()
        >>> len(cs)
        0
        """
        self.thisptr.clear()

    def __reduce__(self):
        """
        Pickle object.

        Examples:

        >>> from ppl import Variable, Congruence_System
        >>> from pickle import loads, dumps
        >>> x = Variable(0)
        >>> y = Variable(1)
        >>> cs = Congruence_System([(3*x+2*y+1 == 3) % 15, (x+2*y + 7 == 0) % 22])
        >>> cs
        Congruence_System {3*x0+2*x1+13==0 (mod 15), x0+2*x1+7==0 (mod 22)}
        >>> loads(dumps(cs))
        Congruence_System {3*x0+2*x1+13==0 (mod 15), x0+2*x1+7==0 (mod 22)}
        """
        return (Congruence_System, (tuple(self),))


cdef _wrap_Congruence(PPL_Congruence congruence):
    cdef Congruence c = Congruence.__new__(Congruence)
    c.thisptr = new PPL_Congruence(congruence)
    return c


def congruence(le, m):
    return Congruence(le, m)
