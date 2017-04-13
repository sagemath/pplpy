# distutils: language = c
# distutils: libraries = gmp

cdef extern from "Python.h":
    ctypedef struct PyObject:
        pass

cdef extern from "gmp.h":
    void mpz_set_si(mpz_t rop, signed long int op)
    void mpz_set(mpz_t rop, mpz_t op)
    
cdef extern from "gmpy2.h":
    cdef (PyObject *)GMPy_MPZ_New(void *)
    cdef int import_gmpy2()
    
cdef extern from "gmpy2_mpz.h":
    cdef mpz_t MPZ(PyObject *)

cdef extern from "gmpy2_types.h":
    ctypedef struct MPZ_Object:
        pass
    ctypedef struct CTXT_Object:
        pass
    
cdef extern from "gmpy2_convert.h":
    # mpz_set_PyStr converts a Python "string" into a mpz_t structure. It accepts
    # a sequence of bytes (i.e. str in Python 2, bytes in Python 3) or a Unicode
    # string (i.e. unicode in Python 3, str in Python 3). Returns -1 on error,
    # 1 if successful.
    cdef int mpz_set_PyStr(mpz_ptr z, PyObject *s, int base);
    
cdef extern from "gmpy2_convert_gmp.h":
    cdef (PyObject *) GMPy_PyLong_From_MPZ(MPZ_Object *obj, CTXT_Object *context)

import_gmpy2()

cdef get_gmpy_mpz(mpz_srcptr z):
    cdef PyObject * res = GMPy_MPZ_New(NULL)
    #cdef PyObject * res = GMPy_MPZ_New(NULL)[0]
    mpz_set(MPZ(res), z)
    # print('bob')
    # print(res[0])
    return <object>res