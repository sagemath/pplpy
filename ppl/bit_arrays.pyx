# distutils: language = c++
# distutils: libraries = gmp ppl m
#*****************************************************************************
#       Copyright (C) 2018 Vincent Delecroix <vincent.delecroix@labri.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from libc.limits cimport ULONG_MAX

from .ppl_decl cimport PPL_Bit_Row

cdef class Bit_Row(object):
    r"""
    A bit row

    This can be considered as a subset of the non-negative integers (with
    upper bound the maximum value ``ULONG_MAX`` that fits in an
    ``unsigned long``).
    """
    def __cinit__(self):
        self.thisptr = new PPL_Bit_Row()

    def __dealloc__(self):
        del self.thisptr

    def __repr__(self):
        r"""
        Examples:

        >>> from ppl import Bit_Row
        >>> r = Bit_Row()
        >>> r
        {}
        >>> r.set(2); r.set(12)
        >>> r
        {2, 12}
        """
        cdef list l = []
        cdef unsigned long k = self.thisptr.first()
        while k != ULONG_MAX:
            l.append(str(k))
            k = self.thisptr.next(k)
        return "{" + ", ".join(l) + "}"

    def __hash__(self):
        r"""
        Examples:

        >>> from ppl import Bit_Row
        >>> r = Bit_Row()
        >>> hash(r)
        Traceback (most recent call last):
        ...
        TypeError: Bit_Row unhashable
        """
        raise TypeError("Bit_Row unhashable")

    def set(self, unsigned long k):
        r"""
        Set the bit in position ``k``

        Examples:

        >>> from ppl import Bit_Row
        >>> r = Bit_Row()
        >>> r.set(2); r.set(34); r.set(12)
        >>> r
        {2, 12, 34}
        >>> r.set(22)
        >>> r
        {2, 12, 22, 34}
        """
        self.thisptr.set(k)

    def set_until(self, unsigned long k):
        r"""
        Set bits up to position ``k`` (excluded)

        Examples:

        >>> from ppl import Bit_Row
        >>> r = Bit_Row()
        >>> r.set_until(5)
        >>> r
        {0, 1, 2, 3, 4}
        """
        self.thisptr.set_until(k)

    def clear_from(self, unsigned long k):
        r"""
        Clear bits from position ``k`` (included) onward

        Examples:

        >>> from ppl import Bit_Row
        >>> r = Bit_Row()
        >>> r.set_until(10)
        >>> r.clear_from(5)
        >>> r
        {0, 1, 2, 3, 4}
        """
        self.thisptr.clear_from(k)

    def clear(self):
        r"""
        Clear all the bits of the row

        Examples:

        >>> from ppl import Bit_Row
        >>> r = Bit_Row()
        >>> r.set(1); r.set(3)
        >>> r
        {1, 3}
        >>> r.clear()
        >>> r
        {}
        """
        self.thisptr.clear()

    def first(self):
        r"""
        Return the index of the first set bit or ``-1`` if no bit is set.

        Examples:

        >>> from ppl import Bit_Row
        >>> r = Bit_Row()
        >>> r.first() == -1
        True
        >>> r.set(127)
        >>> r.first() == 127
        True
        >>> r.set(2)
        >>> r.first() == 2
        True
        >>> r.set(253)
        >>> r.first() == 2
        True
        >>> r.clear()
        >>> r.first() == -1
        True
        """
        cdef unsigned long k = self.thisptr.first()
        if k == ULONG_MAX:
            return -1
        else:
            return k

    def last(self):
        r"""
        Return the index of the last set bit or ``-1`` if no bit is set.

        Examples:

        >>> from ppl import Bit_Row
        >>> r = Bit_Row()
        >>> r.last() == -1
        True
        >>> r.set(127)
        >>> r.last() == 127
        True
        >>> r.set(2)
        >>> r.last() == 127
        True
        >>> r.set(253)
        >>> r.last() == 253
        True
        >>> r.clear()
        >>> r.last() == -1
        True
        """
        cdef unsigned long k = self.thisptr.last()
        if k == ULONG_MAX:
            return -1
        else:
            return k

cdef class Bit_Matrix(object):
    pass
