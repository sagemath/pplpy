# distutils: language = c
# distutils: libraries = gmp

cdef extern from "gmp.h":
    void mpz_set(mpz_t rop, mpz_t op)

import_gmpy2()

cdef GMPy_MPZ_From_mpz(mpz_srcptr z):
    cdef PyObject * res = GMPy_MPZ_New(NULL)
    mpz_set(MPZ(res), z)
    return <object>res

cdef GMPy_MPQ_From_mpq(mpq_srcptr q):
    pass