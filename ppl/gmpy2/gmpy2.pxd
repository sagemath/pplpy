cdef extern from "gmp.h":
    ctypedef struct __mpz_struct:
        pass
    ctypedef __mpz_struct mpz_t[1]
    ctypedef __mpz_struct *mpz_ptr
    ctypedef const __mpz_struct *mpz_srcptr

cdef mpz_get_pylong(mpz_srcptr z)
cdef mpz_get_pyintlong(mpz_srcptr z)