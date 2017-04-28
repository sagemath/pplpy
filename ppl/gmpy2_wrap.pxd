from cpython.object cimport PyObject

cdef extern from "gmp.h":
    # gmp's integer numbers structs
    ctypedef struct __mpz_struct:
        pass
    ctypedef __mpz_struct mpz_t[1]
    ctypedef __mpz_struct *mpz_ptr
    ctypedef const __mpz_struct *mpz_srcptr

    # gmp's rationnal numbers structs
    ctypedef struct __mpq_struct:
        pass
    ctypedef __mpq_struct mpq_t[1]
    ctypedef __mpq_struct *mpq_ptr
    ctypedef const __mpq_struct *mpq_srcptr

    void mpz_init (mpz_t integer)

cdef GMPy_MPZ_From_mpz(mpz_srcptr z)
cdef GMPy_MPQ_From_mpq(mpq_srcptr q)
cdef GMPy_MPQ_From_mpz(mpz_srcptr num, mpz_srcptr den)
#TODO use GMPy_MPZ_From_PyIntOrLong gmpy2_convert_gmp.h instead of mpz_set_pylong
cdef int mpz_set_pylong(mpz_ptr z, L) except -1
