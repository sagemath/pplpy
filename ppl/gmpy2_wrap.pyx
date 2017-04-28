# distutils: language = c
# distutils: libraries = gmp

from cpython.object cimport Py_SIZE

cdef extern from "gmp.h":
    void mpz_set(mpz_t rop, mpz_t op)
    void mpz_neg (mpz_t rop, mpz_t op)
    void mpz_import (mpz_t rop, size_t count, int order, int size, int endian, size_t nails, void *op)
    void mpq_set(mpq_ptr rop, mpq_srcptr op)
    void mpq_set_num (mpq_t rational, mpz_t numerator)
    void mpq_set_den (mpq_t rational, mpz_t denominator)

cdef extern from "gmpy2.h":
    cdef (PyObject *)GMPy_MPZ_New(void *)
    cdef (PyObject *)GMPy_MPQ_New(void *)
    cdef int import_gmpy2()
    cdef mpz_t MPZ(PyObject *)
    cdef mpq_t MPQ(PyObject *)

cdef extern from "longintrepr.h":
    cdef _PyLong_New(Py_ssize_t s)
    cdef long PyLong_SHIFT
    ctypedef unsigned int digit
    ctypedef struct PyLongObject:
        digit* ob_digit

import_gmpy2()

# Unused bits in every PyLong digit
cdef size_t PyLong_nails = 8*sizeof(digit) - PyLong_SHIFT

cdef GMPy_MPZ_From_mpz(mpz_srcptr z):
    """
    Convert from 'gmp mpz' to 'gmpy2 mpz'
    """
    cdef PyObject * res = GMPy_MPZ_New(NULL)
    mpz_set(MPZ(res), z)
    return <object>res

cdef GMPy_MPQ_From_mpq(mpq_srcptr q):
    """
    Convert from 'gmp mpq' to 'gmpy2 mpq'
    """
    cdef PyObject * res = GMPy_MPQ_New(NULL)
    mpq_set(MPQ(res), q)
    return <object>res

cdef GMPy_MPQ_From_mpz(mpz_srcptr numerator, mpz_srcptr denominator):
    """
    Build a 'gmpy2 mpq' from 'gmp mpz' numerator and denominator
    """
    cdef PyObject * res = GMPy_MPQ_New(NULL)
    mpq_set_num(MPQ(res), numerator)
    mpq_set_den(MPQ(res), denominator)
    return <object>res

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
