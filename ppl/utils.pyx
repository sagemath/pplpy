# distutils: language = c
# distutils: libraries = gmp

from cpython.object cimport Py_SIZE

cdef extern from "gmp.h":
    void mpz_set(mpz_t rop, mpz_t op)
    void mpz_neg (mpz_t rop, mpz_t op)
    void mpz_import (mpz_t rop, size_t count, int order, int size, int endian, size_t nails, void *op)

cdef extern from "longintrepr.h":
    cdef _PyLong_New(Py_ssize_t s)
    cdef long PyLong_SHIFT
    ctypedef unsigned int digit
    ctypedef struct PyLongObject:
        digit* ob_digit

# Unused bits in every PyLong digit
cdef size_t PyLong_nails = 8*sizeof(digit) - PyLong_SHIFT

cdef int mpz_set_pylong(mpz_ptr z, L) except -1:
    """
    Convert a Python ``long`` `L` to an ``mpz``.
    """
    cdef Py_ssize_t pylong_size = Py_SIZE(L)
    if pylong_size < 0:
        pylong_size = -pylong_size
    mpz_import(z, pylong_size, -1, sizeof(digit), 0, PyLong_nails,
            (<PyLongObject*>L).ob_digit)
    if Py_SIZE(L) < 0:
        mpz_neg(z, z)
