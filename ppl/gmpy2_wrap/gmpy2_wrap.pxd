cdef extern from "gmp.h":
    ctypedef struct __mpz_struct:
        pass
    ctypedef __mpz_struct mpz_t[1]
    ctypedef __mpz_struct *mpz_ptr
    ctypedef const __mpz_struct *mpz_srcptr

cdef get_gmpy_mpz(mpz_srcptr z)