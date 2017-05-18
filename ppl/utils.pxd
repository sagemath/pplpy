cdef extern from "gmp.h":
    # gmp integer
    ctypedef struct __mpz_struct:
        pass
    ctypedef __mpz_struct mpz_t[1]
    ctypedef __mpz_struct *mpz_ptr
    ctypedef const __mpz_struct *mpz_srcptr

#TODO: make this available from gmpy2
cdef int mpz_set_pylong(mpz_ptr z, L) except -1
